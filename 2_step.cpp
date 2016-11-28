/* Author: Jingyu LI
 * Created at: 23 Aug, 2016
 * Last Modified at: 12:00 28 Nov, 2016
 * Usage: second step for Dan's method, to extract neighbors of a pair of repeat
 * Input: preName.txt and self_compare.sh file
 * output: 2_step.sh
 * Parameters: 
 * 1. Generated shell file in step 1
 * 2. Genome set directory
 * 3. The output shell file name
 * Usage: g++ 2_step.cpp -o 2_step
 * 	  ./2_step 1_step_shell_file Genome_set_directory 2_step_shell_file_name > log_file_name
 * 	  (e.g.: ./2_step 1_step.sh ./Genomes 2_step.sh > 2_step.log)
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

int main(int argc, char * argv[]){
	ifstream fin1;//preName file
	ifstream fin2;
	ifstream fin3;//blast txt file
	ofstream fout;
	fin1.open("./faNamesPre.txt");
	if(!fin1.is_open()){
		cout<<"cannot open pre names file."<<endl;
		return -1;
	}
	fin2.open(argv[1]);
	if(!fin2.is_open()){
		cout<<"cannot open shell file."<<endl;
		return -1;
	}
	fout.open(argv[3]);
	string faNames[25];
	for(int i=0; i<25; i++){
		getline(fin1, faNames[i]);
	}
	fin1.close();
	string str,sortedstr;
	string folder, sorted_file;
	string strainID,label;
	long pos1, pos3, qstart, qend, sstart, send;
	string cpos1, cpos3, cqstart, cqend, csstart, csend;
	stringstream ss;
	int tag;
	cout<<"start process..."<<endl;
	//getline(fin2, str);//mkdir compare_result
	while(!fin2.eof()){
		getline(fin2, str);
		if(str[0]=='c'&&str[1]=='d'&&str[3]=='c'){// "cd compare_result/......"
			folder=str.substr(str.find_first_of(" ")+1);//compare_result/66/66_1
			label=folder.substr(folder.find_first_of("/")+1);//66/66_1
			label=label.substr(0,label.find_first_of("/"));//66
			fout<<str<<"\n";
			getline(fin2, str);//sh .......
			str=str.substr(str.find_first_of(" ")+1);
			str=str.substr(str.find_first_of(" ")+1);
			strainID=str.substr(0,str.find_first_of(" "));
			str=str.substr(str.find_first_of(" ")+1);
			cpos1=str.substr(0,str.find_first_of(" "));
			ss << cpos1;
			ss >> pos1;
			ss.clear();
			str=str.substr(str.find_first_of(" ")+1);
			cpos3=str.substr(0,str.find_first_of(" "));
			ss << cpos3;
			ss >> pos3;
			ss.clear();
			getline(fin2, str);//second sh ...... line
			getline(fin2, str);//bl2seq line
			getline(fin2, str);//sort result line
			sorted_file = str.substr(str.find_last_of(" ")+1);
			//cout<<sorted_file<<"  ";
			sorted_file=folder+"/"+sorted_file;
			fin3.open(sorted_file.c_str());
			if(!fin3.is_open()){
				cout<<"cannot open sorted file."<<sorted_file<<endl;
				//return -1;
			}
			tag=0;//if tag == 0, delete this file
			getline(fin3, sortedstr);
			getline(fin3, sortedstr);
			getline(fin3, sortedstr);
			while(!fin3.eof()){
				getline(fin3, sortedstr);
				if(sortedstr=="") continue;
				for (int i=0;i<6;i++)
					sortedstr=sortedstr.substr(sortedstr.find_first_of("\t")+1);
				//start with q.start
				//cout<<sortedstr<<endl;
				cqstart=sortedstr.substr(0,sortedstr.find_first_of("\t"));
				ss << cqstart; 
				ss >> qstart;
				ss.clear();
				sortedstr=sortedstr.substr(sortedstr.find_first_of("\t")+1);
				cqend=sortedstr.substr(0,sortedstr.find_first_of("\t"));
				ss << cqend;
				ss >> qend;
				ss.clear();
				sortedstr=sortedstr.substr(sortedstr.find_first_of("\t")+1);
				csstart=sortedstr.substr(0,sortedstr.find_first_of("\t"));
				ss << csstart;
				ss >> sstart;
				ss.clear();
				sortedstr=sortedstr.substr(sortedstr.find_first_of("\t")+1);
				csend=sortedstr.substr(0,sortedstr.find_first_of("\t"));
				ss << csend;
				ss >> send;
				ss.clear();
				if(sstart>300000) continue;
				if(sstart>send) continue;
				cout<<label<<" "<<strainID<<" "<<tag+1<<" "<<pos1+qstart<<" "<<pos1+qend<<" "<<pos3+sstart<<" "<<pos3+send<<endl;
				fout<<"sh ../../../extractSq.sh "<<strainID<<" "<<pos1+qstart<<" "<<pos1+qend<<" ../../../faNamesPre.txt "<<argv[2]<<" > "<<label<<"_"<<strainID<<"_A"<<tag+1<<".fa\n";
				fout<<"sh ../../../extractSq.sh "<<strainID<<" "<<pos1+qstart-10000<<" "<<pos1+qstart<<" ../../../faNamesPre.txt "<<argv[2]<<" > "<<label<<"_"<<strainID<<"_LN"<<tag+1<<".fa\n";
				fout<<"sh ../../../extractSq.sh "<<strainID<<" "<<pos3+send<<" "<<pos3+send+10000<<" ../../../faNamesPre.txt "<<argv[2]<<" > "<<label<<"_"<<strainID<<"_RN"<<tag+1<<".fa\n";
				int strain;
				ss << strainID;
				ss >> strain; 
				ss.clear();
				for(int i=0;i<25;i++){
					if(i == strain-1) continue;
					fout<<"../../../blast-2.2.26/bin/bl2seq -i "<<label<<"_"<<strainID<<"_A"<<tag+1<<".fa -j "<<argv[2]<<faNames[i]<<"* -p blastn -e 1e-10 -D 1 -o "<<label<<"_"<<strainID<<"_A"<<tag+1<<"_on"<<i+1<<".txt\n";
					fout<<"../../../blast-2.2.26/bin/bl2seq -i "<<label<<"_"<<strainID<<"_LN"<<tag+1<<".fa -j "<<argv[2]<<faNames[i]<<"* -p blastn -e 1e-10 -D 1 -o "<<label<<"_"<<strainID<<"_LN"<<tag+1<<"_on"<<i+1<<".txt\n";
					fout<<"../../../blast-2.2.26/bin/bl2seq -i "<<label<<"_"<<strainID<<"_RN"<<tag+1<<".fa -j "<<argv[2]<<faNames[i]<<"* -p blastn -e 1e-10 -D 1 -o "<<label<<"_"<<strainID<<"_RN"<<tag+1<<"_on"<<i+1<<".txt\n";
					fout<<"../../../sort_result "<<label<<"_"<<strainID<<"_A"<<tag+1<<"_on"<<i+1<<".txt s-"<<label<<"_"<<strainID<<"_A"<<tag+1<<"_on"<<i+1<<".txt\n";
					fout<<"rm "<<label<<"_"<<strainID<<"_A"<<tag+1<<"_on"<<i+1<<".txt\n";
					fout<<"../../../sort_result "<<label<<"_"<<strainID<<"_LN"<<tag+1<<"_on"<<i+1<<".txt s-"<<label<<"_"<<strainID<<"_LN"<<tag+1<<"_on"<<i+1<<".txt\n";
					fout<<"rm "<<label<<"_"<<strainID<<"_LN"<<tag+1<<"_on"<<i+1<<".txt\n";
					fout<<"../../../sort_result "<<label<<"_"<<strainID<<"_RN"<<tag+1<<"_on"<<i+1<<".txt s-"<<label<<"_"<<strainID<<"_RN"<<tag+1<<"_on"<<i+1<<".txt\n";
					fout<<"rm "<<label<<"_"<<strainID<<"_RN"<<tag+1<<"_on"<<i+1<<".txt\n";
				}
				tag++;
			}
			fin3.close();
			getline(fin2, str);//cd ../../..
			fout<<str<<"\n";
			if(tag==0){
				fout<<"rm -r "<<folder<<"\n";
			}
			
		}
	}


	return 0;
}
