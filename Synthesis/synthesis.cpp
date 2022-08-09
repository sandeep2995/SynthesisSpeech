#include<iostream>
#include<stdlib.h>
#include<cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h> 
#include <conio.h>
#include <sstream>
#include "config.h" //include configuration file that helps in global settings

using namespace std;
std::ostringstream oss;

short samples[600000], MaxAmp=0; //Maximum number of samples we can handle and Maximum amplitude desired for scaling
long double ThresholdZCR, DCshift=0; //Initialize DCShift and Threshold ZCR
double TotalEnergy=0, ThresholdEnergy, NoiseEnergy=0, TotalZCR=0; //TotalZCR is only for first few InitFrames
long start=0, stop=0, framecount, samplecount=0; //start and end marker for speech signal YES or NO
string line, filename;
unsigned short q; //number of coefficients (ai's) that need to be found
long double *R, *c, *E, *k, *alpha;

ifstream InpSpeech;
ofstream ScaledInpSpeech,out,cep;

void durbin()
{
	E = new long double[p+1];
	k = new long double[p+1];
	alpha = new long double[p+1]; // for new alpha values
	long double* alp = new long double[p+1]; // for old values
	if(R[0]==0)
	{
		cout << "\nEnergy should not be ZERO\n";
		return;
	}
	E[0]=R[0]; //Initialize E[0] with energy
	for(int i=0;i<=p;i++)
	{
		alpha[i]=0; //initialize all the initial coefficients to 0
		alp[i]=0;
	}
	for(int i=1;i<=p;i++)
	{
		if(i==1)
			k[1]=R[1]/R[0]; //special case when i=1
		else //find k(i) for all other remaining values
		{
			long double sum=0;
			for(int j=1;j<=i-1;j++)
			{
				alp[j]=alpha[j]; // store the old values in alp array
				sum+=alp[j]*R[i-j];
			}
			k[i]=(R[i]-sum)/E[i-1];
		}
		alpha[i]=k[i]; //assign values to coefficients
		for(int j=1;j<=i-1;j++)
			alpha[j]=alp[j]-k[i]*alp[i-j]; //update coefficients from previously obtained values
		E[i]=(1-k[i]*k[i])*E[i-1]; //update E(i)
	}
	cout << "\nLPC Coefficient values\n";
	for(int i=1;i<=p;i++)
		cout << "a" << i << " = " << alpha[i] << "\n";
}

void cepstralcoefficients()
{
	unsigned int m,j;
	c = new long double[q+1];

	c[0]= logl(R[0]); //initial cepstral coefficent computation from energy
	for(m=1;m<=p;m++)
	{
		long double sum=0;
		for(j=1;j<=m-1;j++) //calculate the sum from older cepstral coefficents to compute new cepstral coefficients
			sum+=j*c[j]*alpha[m-j]/m;
		c[m]=alpha[m]+sum; //new cepstral coefficients
	}
	if(m>p) // for our assignment this never get executed as we assume q=p
	{
		for(;m<=q;m++)
		{
			long double sum=0;
			for(j=m-p;j<=m-1;j++)
				sum+=j*c[j]*alpha[m-j]/m;
			c[m]=sum;
		}
	}
	cout << "\nCepstral Coefficient values\n";//write the resuls to file
	for(int i=0;i<=q;i++)
	{
		cep << c[i] << " ";
	}
	cep << endl;
}

void Coefficients()
{

	q=p; // for our assignment we always take q=p

	durbin(); //call this function to evaluate durbin algorithm
	cepstralcoefficients(); //call method to calculate cepstral coefficients

}

void AutoCorrelation()
{
	long N=framesize, ind=0,K,t;
	R = new long double[p+1]; //create array of size p+1 for p coefficients
	long double *x=new long double[N]; // to store all the samples of a frame
	InpSpeech.open("ScaledInpSpeech.txt", ios::in); // open ScaledInpSpeech.txt to read the scaled samples
	if (!InpSpeech)
	{
		cout << "\n**File failed to open**\n\n";
		InpSpeech.clear();
	}

	if (InpSpeech.is_open()) //if file is open then only execute this block
	{
		cep.open("cep.txt");
		if (!cep) //if file is not present pop the error
		{
			cout << "\n**File failed to open**\n\n";
			cep.close();
			return;
		}
		for(K=0;K<framecount;K++) // iterate through all frames
		{

			if(K < start || K > stop) //ignore the non speech frames
				for(int t=0;t<framesize;t++)
				{
					getline (InpSpeech,line); //read a line at a time
				}
			else //process the speech frames
			{
				for(t=0;t<framesize;t++) // read all the sample values into frames
				{
					getline (InpSpeech,line); //read a line at a time
					x[t]=atof(line.c_str());
				}
				for(ind=0;ind<=p;ind++)  //generate p+1 R values to compute cepstral coefficients
				{
					double square=0;
					for(int t=0;t<framesize-ind;t++)
					{
						square+=x[t]*x[t+ind];
					}
					R[ind]=square/N;
				}
				cout << "frame number " << K <<endl;
				Coefficients(); //call this function to compute the coefficients			
			}
		}

		InpSpeech.close(); // close the file
		cep.close();
	}
}


int main()
{

	long i,j;

	oss << recmod << " "<<duration <<" "<<inpw << " " << inpt; //define the string manually, because macro cant expand inside string
	system(oss.str().c_str()); //call the system function to read directly from the mic
	filename=inpt; //name of the recorded speech text file

	InpSpeech.open(filename, ios::in); // open the file to read from it
	if (!InpSpeech) //if file is not present pop the error
	{
		cout << "\n**File failed to open**\n\n";
		InpSpeech.clear();
	}

	out.open("out.txt"); // open the file to write the results
			     //count the number of samples and frames in the file
	if (InpSpeech.is_open()) //if file is open then only execute this block
	{
		while ( !InpSpeech.eof() ) //read till the end of file is reached
		{ 
			getline (InpSpeech,line); //read a line at a time
			samplecount+=1;     //increment the sample count by 1
			if(samplecount > IgnoreSamples + 4) //first 4 lines in text file indiactes type of encoding of the speech,
				{
					samples[samplecount - (IgnoreSamples + 5)] = (short)atoi(line.c_str()); //5=4+1 4-->indicates encoding, 1-->array index starts with 0
					DCshift+=samples[samplecount - (IgnoreSamples + 5)];
					if(abs(samples[samplecount - (IgnoreSamples + 5)])>MaxAmp)
						MaxAmp=abs(samples[samplecount - (IgnoreSamples + 5)]); //MaxAmp contains the magnitude od the maximum valued(absolute) sample
				}
		}
		InpSpeech.close();
	}
	samplecount = samplecount - (IgnoreSamples + 4);
	framecount = samplecount/framesize;
	DCshift=DCshift/samplecount;
	out << "Number of Samples = " << samplecount << "\n"; // print the results in out.txt file
	out << "Number of frames = " << framecount << "\n";
	out << "DC shift needed is " << DCshift << "\n";
	out << "Maximum Amplitude = " << MaxAmp << "\n";
	cout << "Number of Samples = " << samplecount << "\n"; // print the results on standard output
	cout << "Number of frames = " << framecount << "\n";
	cout << "DC shift needed is " << DCshift << "\n";
	cout << "Maximum Amplitude = " << MaxAmp << "\n";

	//Store scaled samples in different file
	ScaledInpSpeech.open("ScaledInpSpeech.txt"); // open ScaledInpSpeech.txt to write the scaled sample values
	for(i=0;i<samplecount;i++)
		ScaledInpSpeech << (samples[i] - DCshift)*Ampscale/MaxAmp << "\n"; // writing the scaled samples to the file ScaledInpSpeech.txt
										   //for each sample we are applying DC shift then doing amplitude normalization such that maximum sample magnitude is Ampscale
	ScaledInpSpeech.close();
	//use the scaled samples


	//ZCR and Energy calculation
	long double AvgZCR[400], AvgEnergy[400];
	for(i=0;i<framecount;i++)
	{
		AvgZCR[i]=0; // Initialize Average ZCR for each frame
		AvgEnergy[i]=0; // Initialize Average Energy for each frame
		for(j=0;j<framesize-1;j++)
		{
			if((samples[i*framesize + j]-DCshift)*Ampscale/MaxAmp*((samples[i*framesize + j + 1]-DCshift)*Ampscale*1.0 < 0)) // If two adjacent samples have opposite sign then increment the avergae ZCR by 1
				AvgZCR[i]+=1;
			AvgEnergy[i]+=1.0*(samples[i*framesize + j]-DCshift)*Ampscale/MaxAmp*(samples[i*framesize + j]-DCshift)*Ampscale/MaxAmp; // energy is calculated as square of amplitude of the sample value
		}
		//out << "ZCR in "<<i+1<<"th frame is "<< AvgZCR[i] << "\n";
		AvgZCR[i]/=framesize;
		//out << "Avergae ZCR in\t\t "<<i+1<<"th frame is "<< AvgZCR[i] << "\n";
		AvgEnergy[i]+=1.0*(samples[i*framesize + j]-DCshift)*Ampscale/MaxAmp*(samples[i*framesize + j]-DCshift)*Ampscale/MaxAmp; // average energy for all the frames after scaled sample values
																	 //out << "Energy in\t "<<i+1<<"th frame is "<< AvgEnergy[i] << "\n";
		AvgEnergy[i]/=framesize; // average energy per frame after scaled sample values 
					 //out << "Average Energy in \t\t\t"<<i+1<<"th frame is "<< AvgEnergy[i] << "\n";
		TotalEnergy+=AvgEnergy[i];
	}
	//calculate noise energy to decide threshold for energy
	for(i=0;i<InitFrames;i++)
	{
		TotalZCR+=AvgZCR[i];
		NoiseEnergy+=AvgEnergy[i];
	}
	ThresholdZCR=TotalZCR/InitFrames;
	NoiseEnergy/=InitFrames;
	ThresholdZCR*=0.9;
	//ThresholdEnergy=TotalEnergy/framecount; //threshold energy is the average of all averaged energies
	//ThresholdEnergy*=0.9;
	ThresholdEnergy=NoiseEnergy*10;
	bool flag=false;
	//start and end marker of speech
	for(i=0;i<framecount-3;i++)
	{
		//if(AvgEnergy[i+1]>ThresholdEnergy)
		/*speech Activity Detection with start and end framing*/
		if(AvgZCR[i+1]<ThresholdZCR || AvgEnergy[i+1]>ThresholdEnergy)
		{
			if(flag == false && AvgEnergy[i+2]>ThresholdEnergy && AvgEnergy[i+3]>ThresholdEnergy) // starting frame will be identified when the energy crosses threshold value
			{
				start = i  ; //i th frame is the starting frame marker
				out << "Starting frame is "<< start+1 <<"th frame and starting sample is "<< (start+1)*framesize <<"\n"; // write starting frame and sample number in out.txt file
				out << "Starting time = " << 1.0*(start + 1)*framesize/samplingrate << " seconds\n" ;// writing starting sample time in seconds to out.txt
				cout << "Starting frame is "<< start+1 <<"th frame and starting sample is "<< (start+1)*framesize <<"\n"; // write starting frame and sample number on standard output
				cout << "Starting time = " << 1.0*(start + 1)*framesize/samplingrate << " seconds\n" ; // writing starting sample time in seconds to standard output
				flag = true;
			}
		}
		else if(flag == true && AvgZCR[i] > ThresholdZCR && AvgEnergy[i] < ThresholdEnergy && AvgEnergy[i-1] < ThresholdEnergy && AvgEnergy[i-2] < ThresholdEnergy)// ending frame will be identified when the energy crosses threshold value and AvgZCR is lessthan ThresholdZCR
		{
			stop = i ;
			out << "Ending frame is "<< stop+1 <<"th frame and Ending sample is "<< (stop+1)*framesize<<"\n"; // write ending frame and sample number in out.txt file
			out << "Ending time = " << 1.0*(stop + 1)*framesize/samplingrate << " seconds\n" ; // writing ending sample time in seconds to out.txt
			cout << "Ending frame is "<< stop+1 <<"th frame and Ending sample is "<< (stop+1)*framesize<<"\n"; // write ending frame and sample number on standard output
			cout << "Ending time = " << 1.0*(stop + 1)*framesize/samplingrate << " seconds\n" ; // writing ending sample time in seconds to standard output
			flag = false;
			break;
		}
	}

	if(start==0 || stop==0) // condition to check for speech activity
	{
		cout << "No speech Activity"<<endl;
		system("pause"); //pause to see the output
		return 0;
	}

	double testavgE=0, testavgZ=0;
	for(i=start;i<=stop;i++)   
	{
		testavgE+=AvgEnergy[i]; // sum the avergae energy from start to stop frame to compute over all avegare energy.
		testavgZ+=AvgZCR[i]; // sum the avergae zcr from start to stop frame to compute over all avegare ZCR. 
	}
	testavgE/=(stop-start+1); // average energy per sample
	testavgZ/=(stop-start+1); // average ZCR per sample
	out << "Average Energy of the speech signal is "<< testavgE << "\n";  // write average energy to file out.txt
	out << "Average ZCR of the speech signal is "<< testavgZ*framesize << " \n"; // write average energy to standard output
	cout << "Average Energy of the speech signal is "<< testavgE << "\n"; //write average energy to standard output
	cout << "Average ZCR of the speech signal is "<< testavgZ*framesize << " \n"; // write average ZCR to standard output

	AutoCorrelation();
	//Coefficients(); //call this function to compute the coefficients

	cout <<"open cep.txt to see the cepstral coefficients for each frame" <<endl;

	system("pause");

	return 0;
}
