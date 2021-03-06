/*
    Lie algebra library, last compiled 10/4/2012, written in March 2011.
    
    Provides tools for working with algebras similar to universal enveloping algebras of Lie algebras.
    Specifically, this library assumes it is working with a finite dimensional associative algebra
    with a Lie bracket defined between all basis elements, and provides functions that uses this commutator
    to simplify expressions.
    
    Main Functions:
        LieAlgebra::LieAlgebra(filename)
            Constructs a Lie-like algebra according to description in filename
        LieAlgebra::Simplify(expression)
            "Simplifies" expression, with the property that zero expressions will always be recognized.
        Expression::symmetrize()
            returns the symmetrization of the expression
        LieAlgebra::checkJacobi()
            Checks if described algebra satisfies the Jacobi identity
        LieAlgebra::commutator(Expression 1, Expression 2)    
            Computes the commutator of expression 1 and expression 2
        LieAlgebra::fromString(string)
            Parses a string to extract an expression
        LieAlgebra::SymmetricAlgebra()
            Constructs the symmetric algebra with the same generators as given algebra.
    
    TODO- Strip input of all spaces for constructor of LieAlgebra, to be more user friendly.
        
        
*/
#ifndef __LIEALGEBRA_H__
#define __LIEALGEBRA_H__

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cstdio>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::stringstream;
using std::exception;

class BasisE;
class Term;
class Expression;
class LieAlgebra;

// NOTE: PERMUTATIONS ARE CHECKED USING A HASH FUNCTION. CHECK HASH FUNCTION'S ACCURACY!!!!!!!!!

////////////////////////////////////////////////////////////
// Error classes

class FileNotFound: public exception{
    public:
        virtual const char* what() const throw(){
            return "File not found.";
        }    
};

class FormatError: public exception{
    public:
        virtual const char* what() const throw(){
        return "Incorrect file format.";
        }    
};

class NoSuchBasis: public exception{
    public:
        virtual const char* what() const throw(){
            return "Invalid basis vector";
        }
};

class InvalidCoef: public exception{
    public:
        virtual const char* what() const throw(){
            return "Invalid coefficient";
        }    
};

class InvalidExpression: exception{
    public:
        virtual const char* what() const throw(){
            return "Invalid expression";
        }    
};


///////////////////////////////////////////////////////////////
////// Auxilary Functions
bool isDigit(char k){
     if ((k=='0')||(k=='1')||(k=='2')||(k=='3')||(k=='4')||(k=='5')||(k=='6')||(k=='7')||(k=='8')||(k=='9')) return true;
     return false;
}

bool safe_getline(FILE * fp, char* buf) {
    int c;
    int i = 0;
    while ((c = fgetc(fp)) != EOF && c != '\r' && c != '\n') buf[i++] = c;
    buf[i] = 0;
    if (c == EOF) return false;
    return true;
}


// Gets the coefficient (8a*b*c => 8). term is mutated.
void coef(const string term, double* coef, string* remainder){
    stringstream ss;
    ss.str(term);
    int i = 0;
    if (term[i] == '-'){
        *coef = -1;
        i += 1;
    }
    else {
        if ('0' <= term[i] && term[i] <= '9') *coef = 0;
        else *coef = 1;
    }
    while ('0' <= term[i] && term[i] <= '9') {
        *coef = (*coef) * 10 + (term[i] - '0');
        i++;
    }
    if (term[i] == '.') {
        i++;
        double j = 10;
        while ('0' <= term[i] && term[i] <= '9') {
            *coef += term[i]/j;
            i++;
            j *= 10;
        }
    }
    *remainder = term.substr(i);
}

void split(string s, char c, vector<string> * result) {
    stringstream ss;
    for (int i = 0; i < s.length(); i++) {
        if (s[i] == c) {
            result->push_back(ss.str());
            ss.str("");
        }
        else {
            ss << s[i];
        }
    }
    if (ss.str() != "") result->push_back(ss.str());
}

///////////////////////////////////////////////////////////////
/// Main Code
///////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
////// Basis Elements
class BasisE{
      private:
             int id;
             string symbol;
      public:
             friend class LieAlgebra;
             BasisE(){}
             string toString(){
                    return symbol;
             }                 
             bool operator<(BasisE& a){
                  if (symbol<a.toString()) return true;
                  return false;
             }
             bool operator>(BasisE& a){
                  if (symbol>a.toString()) return true;
                  return false;
             }    
};

////////////////////////////////////////////////////////////////
////// Terms
////////////////////////////////////////////////////////////////

class Term{
    public:
        vector<BasisE> TList;
        double coef;
        friend class LieAlgebra;
        
        Term(){
            coef=0;
        }
        
        Term(BasisE x){
            TList.push_back(x);
            coef=1;
        }

        string toString(){
            string symb;
            vector<BasisE>::iterator it;
            
            if (coef==0) return "0";
            else if (coef==-1) symb="-";
            else if (coef==1) symb="";
            else{
                stringstream sstrm;
                sstrm<<coef;
                symb=sstrm.str();
            }
            for (it=TList.begin();it!=TList.end();it++){
                symb+=it->toString();
                symb=symb+"*";
            }
            string::iterator sit;
            sit=symb.end()-1;
            if (*sit=='*')symb.erase(sit);
            return symb;
        }
        
        bool setCoef(double i){
            coef=i;
        }
        
        double getCoef(){
            return coef;
        }
         
        // Operators
        
        Term operator-(){
            Term temp=*this;
            temp.coef=-coef;
            return temp;
        } 
        
        Term operator*=(const BasisE& rhs){
            TList.push_back(rhs);
            return *this;
        }
        
        Term operator*(const BasisE& rhs){
            Term temp=*this;
            temp*=rhs;
            return temp;
        }
         
        Term operator*=(const Term& rhs){
            vector<BasisE>::const_iterator it;
            for (it=rhs.TList.begin();it!=rhs.TList.end();it++){
                TList.push_back(*it);
            }
            coef=coef*rhs.coef;
            reduce();
            return *this;
        }
        
        Term operator*(const Term& rhs){
            Term temp=*this;
            temp*=rhs;
            return temp;
        }
        
        Term operator*=(double i){
            coef*=i;
            reduce();
            return *this;
        }
         
        Term operator*(double i){
            Term temp=*this;
            temp*=i;
            return temp;
        }
        
        bool operator==(Term& rhs){
            if (rhs.coef==coef){
                if (*this|rhs) return true;
                return false;
            }
            return false;
        }
         
        bool operator|(Term& rhs){
            vector<BasisE>::iterator it1, it2;
            rhs.reduce();
            reduce();
            it2=rhs.TList.begin();
            if (rhs.TList.size()!=TList.size()) return false;

            for (it1=TList.begin();it1!=TList.end();it1++){
                if ((it2->toString()!=it1->toString()))return false;
                it2++;
            }
            return true;
        }
        

        bool reduce(){
            if (coef==0){
                vector<BasisE>::iterator it;
                it=TList.begin();
                while (it!=TList.end()){        
                    it=TList.erase(it);
                }
                return true;
            }
            return true;
        }

        Term operator=(const Term& rhs){
            TList=rhs.TList;
            coef=rhs.coef;
            return *this;
        }
        
        // reorders the elements of term according to alphabetical order. The algorithm used is an insertion sort.
        void reorder(){
            // a.compare(b) is negative if a<b
            for (int i=0;i<TList.size();i++){
                for (int j=0;j<i;j++){
                    if (TList[i]<TList[j]){
                        BasisE temp=TList[i];
                        for (int k=i;k>j;k--){
                            TList[k]=TList[k-1];
                        }
                        TList[j]=temp;
                        break;
                    }
                }
            }
        }
    
        
        bool operator<(const Term& other){
            return TList.size()<other.TList.size();
        }
        
    private:
        int hash(){
            vector<BasisE>::iterator it;
            int sum=0;
            for (it=TList.begin();it!=TList.end();it++){
                int k=0,h=0;
                const char* cur=it->toString().c_str();
                while (*cur!=0){
                    h+=*(cur++)<<(k++);
                }
                sum+=h*h;
            }
            return sum;
        }
};



class Expression{
    private:
        // Convert exp to standard form (a+-b=a-b)
        string toStd(string exp){
            string::iterator it;
            for(it=exp.begin();it!=exp.end();it++){
                if ((*it)=='+'){
                    if (*(it+1)!='-') continue;
                    it=exp.erase(it);
                }
            }
            return exp;
        }
        
        // Term a*b*c -> 1/6(a*b*c+a*c*b+b*a*c+b*c*a+c*a*b+c*b*a)
        Expression symmetrize(Term a){
            // First we sort the elements of a
            a.reorder();
            Expression ans;
            int n=a.TList.size();
            ans=ans+a;
            while (1){
                int k,l;
                bool to_break=true;
                for (k=n-2;k>=0;k--){
                    if (a.TList[k]<a.TList[k+1]){
                        to_break=false;
                        break;
                    }
                }
                if (to_break) break;
                for (l=n-1;l>=0;l--){
                    if (a.TList[k]<a.TList[l]) break;
                }
                BasisE temp=a.TList[k];
                a.TList[k]=a.TList[l];
                a.TList[l]=temp;
                for (int j=0;j<(n-k-1)/2;j++){
                    int first=k+1+j;
                    int second=n-1-j;
                    temp=a.TList[first];
                    a.TList[first]=a.TList[second];
                    a.TList[second]=temp;
                }
                ans=ans+a;
            }
            double inversecoef=(double)ans.TList.size();
            ans=ans*(1.0/inversecoef);
            return ans;
        }
    public:
        vector<Term> TList;
        friend class LieAlgebra;
        Expression(){}
      
        Expression(Term x){
            TList.push_back(x);
        }
        Expression(BasisE x){
            TList.push_back(Term(x));
        }
        
        Term getTerm(int i){
            return TList[i];
        }
        
        string toString(){
            string symb;
            vector<Term>::iterator it;
            if (TList.empty()) return "0";
            for (it=TList.begin();it!=TList.end();it++){
                symb=symb+it->toString()+"+";
            }  
            string::iterator sit;
            sit=symb.end()-1;
            symb.erase(sit);
            return toStd(symb);
        }
        // Operators:
        // Addition
        Expression operator+=(const Expression &rhs){
            vector<Term>::const_iterator it;
            for (it=rhs.TList.begin();it!=rhs.TList.end();it++){
                TList.push_back(*it);
            }
            return *this;
        }
        
        Expression operator+(const Expression& rhs){
            Expression temp=*this;
            temp+=rhs;
            return temp;
        }
        
        Expression operator+=(const Term& rhs){
            TList.push_back(rhs);
            return *this;
        }
        
        Expression operator+(const Term& rhs){
            Expression temp=*this;
            temp+=rhs;
            return temp;
        }
        
        // Multiplication
        
        Expression operator*=(const Term& rhs){
            vector<Term>::iterator it;
            for (it=TList.begin();it!=TList.end();it++){
                *it=(*it)*(rhs);
            }
            return *this;
        }
        
        Expression operator*(const Term &rhs){
            Expression temp=*this;
            temp*=rhs;
            return temp;
        }
        
        Expression operator*=(double rhs){
            vector<Term>::iterator it;
            for (it=TList.begin();it!=TList.end();it++){
                *it=(*it)*(rhs);
            }
            return *this;
        }
        
        Expression operator*(double rhs){
            Expression temp=*this;
            temp*=rhs;
            return temp;
        }
        
        Expression operator*(const Expression& rhs){
            Expression temp;
            Term t1,t2;
            vector<Term>::const_iterator it1, it2;
            for (it1=TList.begin();it1!=TList.end();it1++){
                for (it2=rhs.TList.begin();it2!=rhs.TList.end();it2++){
                    t1=*it1;
                    t2=*it2;
                    temp+=t1*t2;
                }
            }
            temp.eliminate();
            return temp;
        }

        Expression operator=(const Expression& rhs){
            TList=rhs.TList;
            return *this;
        }
           
        Expression operator-(){
            vector<Term>::iterator it;
            Expression ans;
            for (it=TList.begin();it!=TList.end();it++){
                ans.TList.push_back(-(*it));
            }
            return ans;
        }
        
        Expression operator-(Expression rhs){
            return *this+(-rhs);
        }
        
        Expression operator-(Term rhs){
            return *this+(-rhs);
        }
        
        // Simplification
        
        bool isZero(){
            if (TList.size()==0) return true;
        }
        
        void eliminate(){
            vector<Term>::iterator it1,it2;
            for (it1=TList.begin();it1!=TList.end();it1++){
                for (it2=it1+1;it2!=TList.end();it2++){
                    if (*it1|*it2){
                        it1->setCoef(it1->getCoef()+it2->getCoef());
                        it2=TList.erase(it2);
                        it2--;
                    }
                }
                if (abs(it1->getCoef())<=0.00000001){
                    it1=TList.erase(it1);
                    it1--;
                    // break;
                }
            }                    
        }
        
        // Symmetrization
        
        Expression symmetrize(){
            Expression ans;
            for (int i=0;i<TList.size();i++){
                ans=ans+symmetrize(TList[i]);
            }
            return ans;
        }
};

const Expression ZERO=Expression();
          
class LieAlgebra{
    private:
        BasisE *basis;
        string *names;
        Expression *ctable;
        int size;
        
        Expression& getR(int i,int j){// Gets commutator of basis
            return *(ctable+size*i+j);
        }
        
        bool setR(Expression e, int i, int j){// Sets commutator of basis
            *(ctable+size*i+j)=e;
            return true;
        }
        
        string toComp(string exp){
            string::iterator it;
            for(it=exp.begin();it!=exp.end();it++){
                if ((*it)=='-'){
                    if (it==exp.begin()) continue;
                    if (*(it-1)=='+') continue;
                    it=exp.insert(it,'+');
                }
            }
            return exp;       
        }
        

    public:
        friend class BasisE;
    
        // Constructors and Destructors
        LieAlgebra(){}
        
        ~LieAlgebra(){
            delete[] basis;
            delete[] names;
            delete[] ctable;
        }
        
        // Reads Lie algebra description from file
        LieAlgebra(string filen){
            FILE * file = fopen(filen.c_str(), "r");
                
            if(!file){
                throw FileNotFound();
            }
            
            try{
                    
                char cpos[500];
                if (!safe_getline(file, cpos)) throw FormatError();
                size=atoi(cpos); // gets number of basis elements
                basis=new BasisE[size];
                names=new string[size];
                ctable=new Expression[size*size]; // need change
                    
                for (int i=0;i<size;i++){
                    if (!safe_getline(file, cpos)) throw FormatError();
                    if (strlen(cpos) == 0) {
                        i--;
                        continue;
                    }
                    *(names+i)=string(cpos);
                    (basis+i)->id=i;
                    (basis+i)->symbol=names[i];
                }

                for (int i=1;i<size;i++)
                    for (int j=0;j<i;j++)
                        setR(fromString("0"),j,i);

                while (safe_getline(file, cpos)){
                    char arg1[500], arg2[500], out[500];
                    int i;
                    char * p;
                    p = cpos;

                    // Read in line of the form \b*[\b*\w+\b*,\b*\w+\b*]\b*=\b*.*
                    if (strlen(p) == 0) continue;
                    while (*(p++) == ' ');
                    p--;
                    if (*(p++) != '[') throw FormatError();
                    while (*(p++) == ' ');
                    p--;
                    i = 0;
                    while (*p != 0 && *p != ',' &&  *p != ' ') {
                        arg1[i++] = *(p++);
                    }
                    while (*(p++) == ' ');
                    p--;
                    arg1[i] = 0;
                    if (*p != ',') throw FormatError();
                    p++;

                    while (*(p++) == ' ');
                    p--;
                    i = 0;
                    while (*p != 0 && *p != ']' &&  *p != ' ') {
                        arg2[i++] = *(p++);
                    }
                    while (*(p++) == ' ');
                    p--;
                    arg2[i] = 0;
                    if (*p != ']') throw FormatError();
                    p++;
                    while (*(p++) == ' ');
                    p--;
                    if (*p != '=') throw FormatError();
                    p++;
                    strcpy(out, p);

                    // Now we have arguments and output
                    int i1=getBasisRef(arg1);
                    int i2=getBasisRef(arg2);
                    if (i1>i2) setR(-fromString(out),min(i1,i2),max(i1,i2));
                    else if (i1<i2) setR(fromString(out),min(i1,i2),max(i1,i2));
                    
                }

                fclose(file);
            }
            catch(...){
                fclose(file);
                throw FormatError();
            }     
        }
        
        // retrieval functions for basis elements
        BasisE& getBasisE(string name){
            return *(basis+getBasisRef(name));
        }
             
        int getBasisRef(string name){
            for (int i=0;i<size;i++){
                if (*(names+i)==name) return i;
            }
            throw NoSuchBasis();
        }
        
        // returns Lie algebra with same basis but in which all commutators are zero
        LieAlgebra SymmetricAlgebra(){
                LieAlgebra g1;
                g1.size=size;
                Expression zero;
                g1.basis=new BasisE[size];
                g1.names=new string[size];
                g1.ctable=new Expression[size*size];// need change
                for (int i=0;i<size;i++){
                    *(g1.names+i)=*(names+i);
                    (g1.basis+i)->id=(basis+i)->id;
                    (g1.basis+i)->symbol=(basis+i)->symbol;
                }
                for (int i=0;i<size;i++){
                    for (int j=i+1;j<size;j++){
                        g1.setR(zero,i,j);
                    }
                }
                return g1;
        }
        
    public:
        // Converts string to expression
        // As this depends on the Lie algebra description, it is a method of the Lie algebra, not expression.
        // For instance, a*b is a valid expression only if the lie algebra description includes a and b.
        Expression fromString(string exp){
            exp=toComp(exp);
            try{
                Term curT;
                Expression curE;
                Expression curParen;
                string curBstr="";
                string curTstr="";
                string::iterator it;
                bool recurse=true;
                double a;

                stringstream transformer;

                // Remove spaces
                for (int i = 0; i < exp.length(); i++) {
                    if (exp[i] != ' ') transformer << exp[i];
                }
                
                exp = transformer.str();

                if (exp == "") return Expression();
                if (exp == "0"){
                    return Expression();
                }

                int i = 0;
                if (exp[0] == '(') {
                    int count = 1;
                    int j = 0;
                    for (j = 1; j < exp.length(); j++) {
                        if (exp[j] == '(') count++;
                        if (exp[j] == ')') count--;
                        if (count == 0) break;
                    }
                    curE = fromString(exp.substr(1, j-1));
                    i = j+1;
                }
                else if (exp[0] != '*' && exp[0] != '+'){
                    int j = 0;
                    int last = 0;
                    double coefficient;
                    string term_str;
                    for (j = 0; j < exp.length(); j++) {
                        if (exp[j] == '(' || exp[j] == '+' || (exp[j] == '-' && j!=0)) break;
                        else if (exp[j] != '*') last = j;
                    }

                    coef(exp.substr(0, last + 1), &coefficient, &term_str);
                    vector<string> basis_elements_str;
                    Term a;
                    a.setCoef(coefficient);
                    split(term_str, '*', &basis_elements_str);
                    for (int k = 0; k < basis_elements_str.size(); k++) {
                        a *= getBasisE(basis_elements_str[k]);
                    }
                    curE = Expression(a);
                    i = last + 1;
                }
                else {
                    throw InvalidExpression();
                }

                if (i == exp.length()) return curE;

                if (exp[i] == '*') return curE * fromString(exp.substr(i+1));
                else if (exp[i] == '+') return curE + fromString(exp.substr(i+1));
                else if (exp[i] == '-') return curE - fromString(exp.substr(i+1));
                else throw InvalidExpression();
            }
            catch(...){
                throw InvalidExpression();
            }
        }        
          
        // commutators  
        Expression commutator(BasisE &x1, BasisE &x2){
            int i1=getBasisRef(x1.symbol);
            int i2=getBasisRef(x2.symbol);
            Expression ans;
            if (i1==i2) return ans;
            if (i1>i2) return -commutator(x2,x1);
            return getR(min(i1,i2),max(i1,i2));
        }
        
        Expression commutator(Expression x, Expression y){
            try{
                return Simplify(x*y-y*x);
            }
            catch(...){
                throw InvalidExpression();
            }
        }
        
        // Poisson bracket ({a*b,c}=a*{b,c}+{a,c}*b, and {a,b}=[a,b] for basis elements)
        Expression poisson(Term x, BasisE& y){
            vector<BasisE>::iterator it;
            Expression ans;
            for (it=x.TList.begin(); it!=x.TList.end();it++){
                Term temp;
                vector<BasisE>::iterator it1;
                for (it1=x.TList.begin(); it1!=x.TList.end();it1++){
                    if (it1==it) continue;
                    temp.TList.push_back(*it1);
                }
                temp.coef=x.coef;
                ans=ans+Expression(temp)*commutator(*it,y);
            }
            return ans;
        }
        
        Expression poisson(Expression x, BasisE& y){
            vector<Term>::iterator it;
            Expression ans;
            for (it=x.TList.begin();it!=x.TList.end();it++){
                ans+= poisson(*it,y);
            }
            return ans;
        }
            
        
        // Checks if the Jacobi identity is satisified.
        bool checkJacobi(){
            for (int i=0; i<size; i++){
                for (int j=i+1; j<size; j++){
                    for (int k=j+1; k<size; k++){
                        Expression x1=Expression(*(basis+i));
                        Expression x2=Expression(*(basis+j));
                        Expression x3=Expression(*(basis+k));
                        Expression check=commutator(commutator(x1,x2),x3)+commutator(commutator(x2,x3),x1)+commutator(commutator(x3,x1),x2);
                        check=Simplify(check);
                        if (!check.isZero()) return false;
                    }
                }
            }
            return true;
        }
        
        
        // returns the side effect of flipping i-th and j-th basis element in given term. The term itself is not included
        Expression flipwc(Term a, int i, int j){
            int index;
            int temp=i;
            i=min(i,j);
            j=max(temp,j);
            Expression ans;
            Term cur=a;
     
            for (index=0;index<j-i;index++){
                Term h1, h2;
                h1.setCoef(a.getCoef());
                h2.setCoef(1);
                vector<BasisE>::iterator it;
                for (it=cur.TList.begin();it!=(cur.TList.begin()+i+index);it++){
                    h1*=(*it); // Grabs first half
                }
                for (it=cur.TList.begin()+i+index+2;it!=cur.TList.end();it++){
                    h2*=(*it); // Grabs second half
                }
                ans+=(Expression(h1)*commutator(*(cur.TList.begin()+i+index),*(cur.TList.begin()+i+index+1)))*Expression(h2);
                cur=vflip(cur,i+index,i+index+1);
            }
            
            for (index=0; index< j-i-1; index++){
                Term h1, h2;
                h1.setCoef(a.getCoef());
                h2.setCoef(1);
                vector<BasisE>::iterator it;
                for (it=cur.TList.begin();it!=(cur.TList.begin()+j-index-2);it++){
                    h1*=(*it);
                }
                for (it=cur.TList.begin()+j-index;it!=cur.TList.end();it++){
                    h2*=(*it);
                }
                ans+=(Expression(h1)*commutator(*(cur.TList.begin()+j-index-2),*(cur.TList.begin()+j-index-1)))*Expression(h2);
                cur=vflip(cur,j-index-2,j-index-1);
            }    
            ans.eliminate();
            return ans;
        }
        
        // flips i-th and j-th entry of term without regarding side-effects. 
        Term vflip(Term a, int i, int j){
            BasisE temp;
            temp=a.TList[i];
            a.TList[i]=a.TList[j];
            a.TList[j]=temp;
            return a;
        }

        // returns an expression equal to the term (as dictated by the lie algebra), with two basis elements of term swapped in position
        Expression flip(Term a, int i, int j){
            int temp=i;
            i=min(i,j);
            j=max(temp,j);
            return Expression(vflip(a,i,j))+flipwc(a,i,j);
        }
        
        // Applies Lie algebra rules to "simplify" expression"
        // The resulting expression may actually be longer than the original,
        // but no two terms in result will have the same multiplicities of basis elements.
        // In particular, a zero expression will always get simplified to 0.
        Expression Simplify(Expression a){
            vector<Term>::iterator Tit1, Tit2;
            
            a.eliminate(); // First get rid of easy stuff
            
            // static int bad=0;
            // static int loops=0;
            try{
            for (Tit1=a.TList.begin();Tit1!=a.TList.end();Tit1++){
                for (Tit2=Tit1+1; Tit2!=a.TList.end();Tit2++){
                    // loops++;
                    if (Tit2->TList.size()!=Tit1->TList.size()) continue; // Easy way out
                    // shouldskip=(Tit2->hash()!=Tit1->hash());
                    if (Tit2->hash()!=Tit1->hash()) continue; // another way out
                    vector<BasisE>::iterator Bit1,Bit2;
                    vector<int> Permutation;
                    bool added_list[Tit2->TList.size()];
                    for (int i=0;i<Tit2->TList.size();i++){
                        added_list[i]=false;
                    }
                    int i1, next;
                    bool isPerm=true;
                    
                    for (Bit1=Tit1->TList.begin();Bit1!=Tit1->TList.end();Bit1++){
                        i1=0;
                        bool found=false;
                        for (Bit2=Tit2->TList.begin();Bit2!=Tit2->TList.end();Bit2++){
                            if (Bit2->toString()==Bit1->toString()){
                                vector<int>::iterator Iit;
                                if (!added_list[i1]){
                                    Permutation.push_back(i1);
                                    found=true;
                                    added_list[i1]=true;
                                    break;
                                }
                            }
                            i1++;
                        }
                        if (found) continue;
                        isPerm=false;
                        break;
                    }
                    if (!isPerm) continue;
            // Now, we know it is permutation
            // What does Permutation tell us? That the ith of the first term gets sent to Permutation[i]th of second
            
                    Expression newTerms;
                    Term curterm=*Tit1;
                    int curperm[Tit1->TList.size()];
                    
                    for (int j=0;j<Tit1->TList.size();j++){
                        curperm[j]=j;
                    }
                    //////////////////////////////////////////////
                    ///// Explanation of the following Loop
                    //////////////////////////////////////////////
                    /*
                        The curperm[i] tracks the original location of the elements of the first term.
                        More exactly, curperm[i] tells the original location of ith term.
                        Say the first term was x1*x2*x3*x4*x5 and the second term was x2*x4*x5*x3*x1
                        We start with the first element of term1, which is curterm.
                        We check where it is supposed to be-It is supposed to be at Permutation[0] of second term
                        Are the 0th and Permutation[0]th elements of the first term equal? 
                        If they are, we do not need to do anything.Otherwise, let us flip the first term.
                        In this case, they are not equal, so we flip first term from 0 to Permutation[curperm[0]]=4. What do we get?
                        We get x5*x2*x3*x4*x1+other terms. Let us ignore the other terms for sake of argument.
                        We set curperm[0]=4
                        and curperm[4]=0
                        We check again. Is the current 0th element equal to Permutation[curperm[0]=4]=2 of current? No. 
                        So we get x3*x2*x5*x4*x1
                        We set curperm[0]=2
                        curperm[2]=4
                        ->0=Permutation[curperm[0]]->curterm[0]=curterm[Permutation[2]=3]? no.
                        curterm->x4*x2*x5*x3*x1
                        curperm[0]=3
                        curperm[3]=2
                        ->curterm[0]=curterm[Permutation[curperm[0]]]? No., ...etc
                    */
                    for (int i=0; i<Tit1->TList.size();i++){
                        for(;;){
                            int t1;
                            next=Permutation[curperm[i]];
                            if (curterm.TList[next].toString()==curterm.TList[i].toString()) break;
                            newTerms=newTerms+flipwc(curterm,min(next,i),max(next,i));
                            curterm=vflip(curterm,min(next,i),max(next,i));
                            t1=curperm[i];
                            curperm[i]=curperm[next];
                            curperm[next]=t1;
                        }
                    }
                    *Tit1=curterm;
                    return Simplify(a+newTerms);
                }
            }
            }
            catch(...){
                return a;
            }
            
            return a;
        }
        // Checks if given expression is central in lie algebra.
        bool isCentral(Expression z){
            string result;
            BasisE j;
            bool answer=true;
            for (int i=0;i<size;i++){
                j=*(basis+i);
                result=commutator(z,Expression(Term(j))).toString();
                cout<<"[element,"<<j.toString()<<"] = "<<result<<endl<<endl;
                if (result!="0") answer=false;
            }
            if (answer) cout<<endl<<z.toString()<<" is IN the center"<<endl;
            else cout<<endl<<z.toString()<<" is NOT in the center"<<endl;
            return answer;
        }
};

#endif
