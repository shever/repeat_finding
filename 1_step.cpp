/* Author: Jingyu Li
 * Created at : 12:00 Thu 18 Aug, 2016
 * Last Modified at: 11:16 Mon 28 Nov, 2016
 * Three Parameters:
 * 1. input maf file
 * 2. Genome data folder directory
 * 3. output shell file (a file name with a .sh format)
 * Usage: g++ 1_step.cpp -o 1_step
 * 	  ./1_step input_maf_file Genome_set_Directory output_shell_file
 * 	  (e.g. ./1_step cordis.maf ./Genomes 1_step.sh)
 * */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char * argv[]){
	ifstream fin1;//preName file
	ifstream fin2;
	ofstream fout;
	fin1.open("./faNamesPre.txt");
	if(!fin1.is_open()){
		cout<<"cannot open pre names file."<<endl;
		return -1;
	}
	fin2.open(argv[1]);
	if(!fin2.is_open()){
		cout<<"cannot open segment file."<<endl;
		return -1;
	}
	fout.open(argv[3]);
	string faNames[25];
	for(int i=0; i<25; i++){
		getline(fin1, faNames[i]);
	}
	fin1.close();
	string str;
	string label;
	stringstream ss;
	string preName;
	long start=0,length=0,strain_len=0;
	int mult,j=0,para1_extract=0;
	char strand;
	string cstart,clength,cstrain_len,strainID;
	string current,preA,target,result,s_result;
	cout<<"start process..."<<endl;
	fout<<"mkdir compare_result\n";
	while(!fin2.eof()){
		getline(fin2, str);
		if(str[0]=='a'){
			str=str.substr(str.find_first_of("=")+1);
			str=str.substr(str.find_first_of("=")+1);
			label=str.substr(0,str.find_first_of(" "));
			fout<<"mkdir compare_result/"<<label<<"\n";
			str=str.substr(str.find_first_of("=")+1);
			ss << str;
			ss >> mult;
			ss.clear();
			j=0;
			for(int i=0;i<mult;i++){
				getline(fin2,str);
				str=str.substr(str.find_first_of(" ")+1);
				preName=str.substr(0, str.find_first_of("."));
				while(preName!=faNames[j] && j<25){
					//fout<<"mkdir compare_result/"<<label<<"/no"<<label<<"_"<<j+1<<"\n";
					j++;
				}
				if(j==25) break;
				str=str.substr(str.find_last_of("\t")+1);//start with start_position
				cstart=str.substr(0,str.find_first_of(" "));
				str=str.substr(str.find_first_of(" ")+1);//start with length
				clength=str.substr(0,str.find_first_of(" "));
				str=str.substr(str.find_first_of(" ")+1);//start with direction
				strand=str[0];
				str=str.substr(str.find_first_of(" ")+1);//start with strain length
				cstrain_len=str.substr(0,str.find_first_of(" "));
				ss << cstart; ss >> start; ss.clear();
				ss << cstrain_len; ss >> strain_len; ss.clear();
				ss << clength;
				ss >> length;
				ss.clear();
				if(strand=='-'){
					start = strain_len-start-length+1;
				}
				j++;
				if(length<2000) continue; 
				ss << j;
				ss >> strainID;
				ss.clear();
				current="compare_result/"+label+"/"+label+"_"+strainID;
				preA=label+strand+"_"+strainID+"_preA.fa";
				target=label+strand+"_"+strainID+"_target.fa";
				result=label+strand+"_"+strainID+".txt";
				s_result="s_"+result;
				fout<<"mkdir "<<current<<"\n";
				fout<<"cd "<<current<<"\n";
				fout<<"sh ../../../extractSq.sh "<<j<<" "<<start-14000<<" "<<start+1000<<" ../../../faNamesPre.txt "<<argv[2]<<" > "<<preA<<"\n";
				fout<<"sh ../../../extractSq.sh "<<j<<" "<<start+1000<<" "<<strain_len<<" ../../../faNamesPre.txt "<<argv[2]<<" > "<<target<<"\n";
				fout<<"../../../blast-2.2.26/bin/bl2seq -i "<<preA<<" -j "<<target<<" -p blastn -e 1e-10 -D 1 -o "<<result<<"\n";
				fout<<"../../../sort_result "<<result<<" "<<s_result<<"\n";
				fout<<"cd ../../..\n";
				
				
			}
		}
	}
	
	return 0;
}
