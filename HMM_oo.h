/*
 *  HMM_oo.h
 *  HMM_object_oriented
 *
 *  Created by Christian on 8/19/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef HMM_OO_H
#define HMM_OO_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "HMM_binner.h"
#include "EmissionProb.h"

using namespace std;

const int HMM_DATA_INT = 1;
const int HMM_DATA_FP = 2;

template <typename T>
class CHMMData
{
private:
	vector<T> m_data;
	vector<short int> m_classifications;
	int m_num_data_points;
public:
	CHMMData();
	void ReadHMMDataFromTxt(string filename,const int data_type);
	void WriteHMMDataToTxt(string filename);			// writes the (classified!!) HMM-Data to a file..
	vector<T>& ReturnRefToData() {return m_data;}
	vector<short int>& RetRefToClassifications(){return m_classifications;}
	int ReturnNumberOfDataPoints(){return m_num_data_points;}
};

template <typename T>
class CHMM
{
private:
	CHMMData<T>* m_p_to_data;
	CBinner<T> m_binner;
	CEmissionPObject<T> m_emission_p_object; 
	
	unsigned short int m_num_bins; // set binning of emission-values
	int m_num_states; // number of states with which to model the data
	int m_order_of_markov_chain; // oder of Markov-Chain, for now 1, might be extended at some point !!!
	vector<double> m_initial_state_probs;
	vector<vector<double> > m_state_transition_probs;
	vector<unsigned short int> m_binned_values;
	vector<vector<double> > m_state_emission_probabilities;
public:
	//CHMM();
	void SetPointerToHMMData(CHMMData<T>* p_to_data){ m_p_to_data=p_to_data;}
	void SetNumberOfStates(int nstates){ m_num_states=nstates;}
	void SetOrderOfMarkovChain(int order){ m_order_of_markov_chain=order;}
	void SetNumberOfBins(int nbins){m_num_bins = nbins;}
	void SetInitialProbs(int nstates,vector<T>& initial_probs);
	void SetTransitionProbs(int nstates,vector<vector<T> >& transition_probs);
	void PrintOutCurrentValues();	// prints out the current values (number of states, initial state probabilities, transition probabilities) to the screen....
	void Iterate();
	void DoBinningOfObservations();	// do the binning of the values...
	void InitializeEmissionProbabilities(const int type, vector<vector<T> >& metadata);
};

template <typename T>
void CHMM<T>::InitializeEmissionProbabilities(const int type, vector<vector<T> >& metadata)
{
	m_emission_p_object.SetNumberOfStates(m_num_states);
	m_emission_p_object.SetType(type);
	if(type==EMISSION_TYPE_GAUSSIAN)
	{
		cout << "Gaussian Emission Probabilities: " << endl;
		m_emission_p_object.SetNumberOfMetadataValues(2);
		m_emission_p_object.SetMetaData(metadata);
		
		// now calculate for each value of the binned values
		// 1. The unbinned value (using m_binner)
		// 2. the emission value for each state (using m_emission_p_object)
		for(int index_state=0;index_state<m_num_states;index_state++)
		{
			vector<double> temp_vec;
			temp_vec.clear();
			double sum=0.0f;
			for(int index_bin=0;index_bin<m_num_bins;index_bin++)
			{
				// get value from index (using m_binner)
				T value = m_binner.ReturnValueFromIndex((unsigned short int)index_bin);
				
				// get emission probability (UNNORMALIZED!!!!)
				// (using m_emission_p_object)
				double prob = m_emission_p_object.ReturnEmissionProbability(index_state,value);
				sum+=prob;
				temp_vec.push_back(prob);
			}
			// normalize and "push_back"
			for(int index_bin=0;index_bin<m_num_bins;index_bin++)
				temp_vec[index_bin]/=sum;
			
			m_state_emission_probabilities.push_back(temp_vec);
		}
	}
}

template <typename T>
void CHMM<T>::DoBinningOfObservations()
{
	// create "binning-object:"
	m_binner.SetNumBins(m_num_bins);
	m_binner.DoBinning(m_p_to_data->ReturnRefToData(),m_binned_values);
	m_binner.PrintOutBinnerParameters();
}

template <typename T>
void CHMM<T>::SetInitialProbs(int nstates,vector<T>& initial_probs)
{
	// clear all existing elements
	m_initial_state_probs.clear();
	for(int i=0;i<nstates;++i)
		m_initial_state_probs.push_back(initial_probs[i]);
}

template <typename T>
void CHMM<T>::SetTransitionProbs(int nstates,vector<vector<T> >& transition_probs)
{
	m_state_transition_probs.clear();
	vector<T> row;
	
	for(int i=0;i<nstates;++i)
	{
		row.clear();
		for(int j=0;j<nstates;++j)
			row.push_back(transition_probs[i][j]);
		
		m_state_transition_probs.push_back(row);
	}
}

template <typename T>
void CHMM<T>::PrintOutCurrentValues()
{
	// print out basic information, like
	// number of states, order of Markov Chain,
	// the initial state probabilities and
	// the state transition probabilities....
	cout << "number of states: " << m_num_states << endl;
	cout << "order of Markov-Chain: " << m_order_of_markov_chain << endl;
	
	cout << "initial state probabilities: " << endl;
	for(int i=0;i<m_num_states;++i)
		cout << "state" << i << ": " << m_initial_state_probs[i] << endl;

	cout << "state transition probabilities: " << endl;
	for(int i=0;i<m_num_states;++i)
	{
		for(int j=0;j<m_num_states;++j)
			cout << m_state_transition_probs[i][j] << "\t";
		
		cout << endl;
	}	 
}

template <typename T> 
CHMMData<T>::CHMMData()
{
	m_num_data_points=0;
}

template <typename T> 
void CHMMData<T>::ReadHMMDataFromTxt(string filename,const int data_type)
{
	string line;
	ifstream data_file(filename.c_str(),ios::in);
	int num_data_points_read=0;
	char* pEnd;
	T temp_value;
	if(data_file.is_open())
	{
		while (getline(data_file,line))
		{
		//	cout << line << endl;
			switch(data_type)
			{
				case HMM_DATA_INT:
					temp_value = strtol(line.c_str(),&pEnd,10);
					break;
				case HMM_DATA_FP: 
					temp_value = strtod(line.c_str(),&pEnd);
					break;
			}
			m_data.push_back(temp_value);
			num_data_points_read++;
		}
	}
	else
		cout << "unable to open file --> check filename" << endl;
	
	m_num_data_points=num_data_points_read;
	cout << num_data_points_read << " data points read from file: " << filename << endl;
	data_file.close();
}

template <typename T> 
void CHMMData<T>::WriteHMMDataToTxt(string filename)
{
	ofstream data_file(filename.c_str(),ios::out);
	if(data_file.is_open())
	{
		for(int i=0;i<m_num_data_points;++i)
		{
			data_file << m_data[i] << "\t" << m_classifications[i] << endl;
		}
	}
	else 
		cout << "unable to open file --> check filename" << endl;
}

#endif