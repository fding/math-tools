// A simple command line interface to the functions of LieAlgebra.h
#include <iostream>
#include "LieAlgebra.h"

using namespace std;


int main(){
    try{
        cout<<"Enter file with algebra description"<<endl;
        string filename;
        cin>>filename;
        LieAlgebra g(filename);
        cout<<endl<<"Loaded algebra in "<<filename;
        
        for (;;){
            cout<<endl<<endl<<"Enter 1 for commutator, 2 for flipping around monomials, 3 to check if expression is central,";
            cout<<"4 to simplify expression, 9 to quit: ";
            int ans;
            stringstream stm;
            string temp;
            cin>>temp;
            stm.str(temp);
            if (!(stm>>ans)) continue;
            if (ans==9) break;
            if (ans==1){
                cout<<endl<<"First expression: ";
                string x1, x2;
                cin>>x1;
                cout<<endl<<"Second expression: ";
                cin>>x2; 
                try{
                    cout<<"["<<x1<<","<<x2<<"] = "<<g.commutator(g.fromString(x1),g.fromString(x2)).toString();
                }
                catch(exception& e){
                    cout<<e.what();
                }
                catch (...){cout<<"Unknown Error";}
            }
            if (ans==2){
                cout<<endl<<"Enter expression: ";
                string expr;
                cin>> expr;
                cout<<endl<<"Which two coordinates?: ";
                cout<<endl;
                cout<<"(if we have expression x1*x2*x3*x4*x5, coordinate of x3 is 3)"<<endl;
                int i1, i2;
                try{
                cin>>i1>>i2;
                i1=i1-1;
                i2=i2-1;
                cout<<expr<<" = "<<g.flip(g.fromString(expr).getTerm(0),i1,i2).toString();
                }
                catch (exception& e){
                      cout<<e.what();
                }
                catch (...){cout<<"Unknown Error";}
                            
            }
            if (ans==3){
               try{
                   string response;
                   cout<<endl<<"Enter expression to check: ";
                   cin>>response;
                   g.isCentral(g.fromString(response));
               }
               catch(exception& e){
                                cout<<e.what();
               }
               catch (...){cout<<"Unknown Error";}
            }
            if (ans==4){
                try{
                    string response;
                    cout<<endl<<"Enter expression to simplify: ";
                    cin>>response;
                    cout<<g.Simplify(g.fromString(response)).toString();
                }
                catch(exception& e){
                                 cout<<e.what();
                }
                catch (...){cout<<"Unknown Error";}
            }
        }
    }
    catch (exception& e){
          cout<<endl<<"Error: "<<e.what()<<endl;
          main();
    }
    
    return 0;
}
