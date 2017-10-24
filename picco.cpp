#ifndef __PICCO__
#define __PICCO__

#include <Rcpp.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include<iostream>
#include<sstream>
#include<string>
#include<vector>
#include<stdlib.h>

// [[Rcpp::depends(BH)]]

using namespace Rcpp;
using namespace std;
using namespace boost; 
typedef boost::adjacency_list<vecS, vecS, undirectedS> mygraph; 

class picco{
	private:
		
		int * value;
		unsigned long int n;
		
	public:
		//Default constructor		
		picco(): value(), n() {}
		
		//Costruttore con valori e lunghezza
		picco(int *count0, unsigned long int x){
			n = x;
			value = new int[n];
			for(int i=0; i<n; i++)
				value[i]=count0[i];
		}
		
		//Copy Constructor
		picco(const picco& p){
			n=p.n;
			value = new int[n];
			
			for(int i=0; i<n; i++)
				value[i]=p.value[i];
		}
		
		//Copy Assignment
		picco& operator= (const picco& p){
		
			if(this!=&p){ //beware of self-assignment
		
				if(n!=p.n){
		
					//deallocazione dell'oggetto
					delete[] value;
			
					//creazione del nuovo oggetto come nel costruttore di copia
					n=p.n;
					value = new int[n];
					for(int i=0; i<n; i++)
						value[i]=p.value[i];

					//warning
					cerr<<"dimensioni del picco non uguali, operazione forzata"<<endl;

				}else{
			
					//modifica dell'oggetto come nel costruttore di copia
					for(int i=0; i<n; i++)
						value[i]=p.value[i];
				}
			
			}
			
			return *this;
		}	
		
		//Distructor
		~picco(){ delete [] value; }
		
		//Peak features
		//Max high
		int maxhigh() const{
		
			int mh=0;
			for(int j=0;j<n;++j)
				if(value[j]>mh)
					mh=value[j];
			
			return mh;
		}
		
		//Area of the peak
		double area() const{
			
			double h=1.0;
			double sum=0.0;
			
			for(int i=0;i<n;++i)
				sum+=value[i]*h;
		
			return sum;
		}
		
		//Maximum matching
		int M() const{
			
			//Min value
			/*int minh=value[0];
			for(int j=0;j<n;++j)
				if(value[j]<minh)
					minh=value[j];
			*/
			
			//Graph underline the peak
			mygraph g;
			
			//long int cursor=value[0]-minh;
			//long int count=cursor+1;
			
			long int cursor=0;
			long int count=1;
			
			mygraph::adjacency_iterator neighbourIt, neighbourEnd;
			
			//
			for(unsigned long int j=1;j<n;++j){
				
				int jump = value[j]-value[j-1];
				
				if(jump>0){
					for(int i=1;i<=jump;++i){
						add_edge(cursor, count, g);
						cursor=count;
						++count;
					}
				}
				else if((jump<0)){
					for(int i=1;i<=-jump;++i){
						tie(neighbourIt, neighbourEnd) = adjacent_vertices(cursor, g);
						for (; neighbourIt != neighbourEnd; ++neighbourIt)
							if(*neighbourIt<cursor)
								cursor=*neighbourIt;
					}
				}
			}
			//	
			
			/*for(unsigned long int j=1;j<n;++j){
				
				int jump = value[j]-value[j-1];
				
				if(jump>0){
					for(int i=1;i<=jump;++i){
						add_edge(cursor, count, g);
						cursor=count;
						++count;
					}
				}
				else if((jump<0)){
					for(int i=1;i<=-jump;++i){
						if(num_vertices(g)!=0){
							int parent=cursor;
					
							tie(neighbourIt, neighbourEnd) = adjacent_vertices(cursor, g);
							for (; neighbourIt != neighbourEnd; ++neighbourIt)
										if(*neighbourIt<cursor)
											parent=*neighbourIt;
					
							if(parent!=cursor){
								cursor=parent;
							}
							else{
								add_edge(cursor, cursor-1, g);
								cursor--;
							}
						}
						else{
							add_edge(cursor,cursor-1,g);
							cursor--;
						}
					}
				}
			}*/
			
			//Print out
			/*mygraph::vertex_iterator vertexIt, vertexEnd;
			
			tie(vertexIt,vertexEnd)=vertices(g);
			
			for(;vertexIt!=vertexEnd;++vertexIt){
				cout << *vertexIt << " is connected with "; 
				tie(neighbourIt, neighbourEnd) = adjacent_vertices(*vertexIt, g); 
				for (; neighbourIt != neighbourEnd; ++neighbourIt) 
    				cout << *neighbourIt << " "; 
    			cout << "\n";
			}*/

			//Maximum matching
			unsigned long int n_vertices=num_vertices(g);
			std::vector<graph_traits<mygraph>::vertex_descriptor> mate(n_vertices);
			
			edmonds_maximum_cardinality_matching(g, &mate[0]);
			
			return matching_size(g, &mate[0]);
		}
		
		unsigned long int getn() const{
			return n;
		}
		
		//Overload <<
		friend ostream& operator<<(ostream& OUT, const picco& p){

			for(unsigned long int i=0;i<p.n;i++)
					OUT<<p.value[i]<<'\t';

			return OUT;
		}
		
		
};
#endif


// [[Rcpp::export]]
int get_M(unsigned long int x,IntegerVector count0){
  std::vector<int> data = Rcpp::as<std::vector<int> >(count0);
	picco peak(&data[0],x);
	int maxmatch=peak.M();
	return maxmatch;
}
