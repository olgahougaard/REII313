#include <QCoreApplication>
#include <iostream>
#include <iomanip>
#include "matrix.h"

using namespace std;

void print_matrix(const Matrix<double> &);
bool Simplex(Matrix<double> &,Matrix<double> &,Matrix<double> &,Matrix<double> &,Matrix<double> &,Matrix<double> &,int,Matrix<double> &,Matrix<double> &);
bool Base(Matrix<double> &,Matrix<double> &,Matrix<double> &,Matrix<double> &,Matrix<double> &,Matrix<double> &,int,Matrix<double> &,Matrix<double> &);
void Branch(Matrix<double> &,Matrix<double> &,Matrix<double> &,Matrix<double> &,Matrix<double> &,int,Matrix<double> &,Matrix<double> &,int);

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    cout<<"*Simplex and Branch & Bound LP Solver*"<<endl;
    bool ob = true;//max=true min=false
    int v, e;
    v=2;
    e=3;

    int iter = 10;//max iterations for simplex
    int range = 1;//range for the branch and bound

    IdentityMatrix<double> B(e,e); //Matrixes
    Matrix<double> xN(v,1);
    Matrix<double> cB(1,e);


    Matrix<double> N_or(e,v);
    Matrix<double> b_or(e,1);
    Matrix<double> cN_or(1,v);

    N_or.put(0,0,0.1); N_or.put(0,1,3);
    N_or.put(1,0,-30); N_or.put(1,1,5);
    N_or.put(2,0,0.6); N_or.put(2,1,1);
    //N_or.put(3,0,-1); N_or.put(3,1,0);


    b_or.put(0,0,5000);
    b_or.put(1,0,5000);
    b_or.put(2,0,4800);
    //b_or.put(3,0,-20);



    cN_or.put(0,0,-60);cN_or.put(0,1,900);


    Matrix<double> N=N_or;
    Matrix<double> b=b_or;
    Matrix<double> cN=cN_or;

    for(int i =0;i<e;i++)
    {
        for(int j=0;j<e;j++)
        {
            if(i==j)
            {
                B.put(i,j,1);
            }
            else
            {
                B.put(i,j,0);
            }
        }
    }

    for(int i=0;i<v;i++)
    {
        xN.put(i,0,0);
    }
    for(int i=0;i<e;i++)
    {
        cB.put(0,i,0);
    }

    Matrix<double> x_B(e,1);//Symbol Matrixes
    Matrix<double> x_N(v,1);
    for(int i=0;i<v;i++)
    {
        x_N.put(i,0,i+1);
    }
    for(int i=0;i<e;i++)
    {
        x_B.put(i,0,v+i+1);
    }


    cout<<endl<<endl;   //Print
    if(ob)
    {
        cout<<"Maximize: ";
    }
    else
    {
        cout<<"Minimize: ";
    }
    for(int i=0;i<v;i++)
    {
        cout<<cN.get(0,i)<<"x"<<i+1;
        if(i!=v-1)
        {
           cout<<" + ";
        }
    }
    cout<<endl<<endl<<"Subject to:"<<endl;
    for(int i=0;i<e;i++)
    {
        for(int j=0;j<v;j++)
        {
            cout<<N.get(i,j)<<"x"<<j+1;
            if(j!=v-1)
            {
               cout<<" + ";
            }
        }
        cout<<" <= "<<b.get(i,0)<<endl;
    }
    cout<<endl<<endl;


    if(Base(B,N,b,xN,cN,cB,0,x_B,x_N))//Find the starting pivot point
    {
        cout<<"Successfully found a starting pivot point";
    }
    else
    {
        cout<<"Could not find a starting pivot point";
    }



    cout<<endl<<"Executing the Simplex Method:"<<endl;
    if(Simplex(B,N,b,xN,cN,cB,iter,x_B,x_N))
    {
        cout<<"Successfully executed the Simplex method"<<endl;
    }
    else
    {
        cout<<"Could not complete the Simplex method"<<endl;
    }


    Matrix<double> xB = B.getInverse()*b -B.getInverse()*N*xN;
    //std::cout << std::setprecision(5) << std::fixed;    //Set decimal place
    cout<<"\nResults: Maximazation Function (Z) = ";
    print_matrix(cB*xB+cN*xN);
    cout<<endl<<"X Varibles values\n";
    cout<<"---------------------------------------------"<<endl;
    /*for(int i=0;i<e+v;i++)
    {
        int j=0;
        bool flag=true;
        for(int j=0;j<x_N.getRows;j++)
        {
            if(i+1==x_N)
            {
            }
        }
        for(int j=0;j<x_B.getRows;j++)
        {
        }
    }*/
    print_matrix(xN);
    print_matrix(x_N);
    cout<<endl;
    print_matrix(xB);
    print_matrix(x_B);
    cout<<"----------------------------------------------"<<endl;

    Matrix<double> x(v,1);
    for(int i=0;i<v;i++)
    {
       if(x_N.get(i,0)<=v)
       {
           x.put(x_N.get(i,0)-1,0,xN.get(i,0));
       }
    }
    for(int i=0;i<e;i++)
    {
       if(x_B.get(i,0)<=v)
       {
           x.put(x_B.get(i,0)-1,0,xB.get(i,0));
       }
    }
    Matrix<double> z(1,1);
    z.put(0,0,0);
    Matrix<double> z2(1,1);
    z2=cB*xB+cN*xN;

    Matrix<double> x2=x;
    Branch(x,x2,N_or,b_or,cN_or,iter,z,z2,range);

    cout<<endl<<"---------------------------------"<<endl;
    print_matrix(x);
    cout<<endl;
    print_matrix(z);
    cout<<endl;

    return a.exec();
}

bool Base(Matrix<double> &B,Matrix<double> &N,Matrix<double> &b,Matrix<double> &xN,Matrix<double> &cN,Matrix<double> &cB,int index,Matrix<double> & x_B,Matrix<double> &x_N)
{
    int e=b.getRows();
    int v=N.getColumns();
    bool flag;
    Matrix<double> xB = B.getInverse()*b -B.getInverse()*N*xN;
    double temp=0;
    Matrix<double> tmp2 = B.getColumn(0);

     flag=true;
         for(int j=0;j<xN.getRows();j++)
         {
             if(xN.get(j,0)<0)
             {
                flag=false;
             }
         }
         for(int j=0;j<xB.getRows();j++)
         {
             if(xB.get(j,0)<0)
             {
                 flag=false;
             }
         }
      if(flag)
      {
          return true;
      }

    for(int i=0;i<e;i++)
    {


         if(index!=v-1)
         {
            if(Base(B,N,b,xN,cN,cB,index+1,x_B,x_N))
            {
                 return true;

            }
         }

         tmp2 = B.getColumn(i);//Swap the columns
         for(int k=0;k<B.getRows();k++)
         {
             B.put(k,i,N.get(k,index));
             N.put(k,index,tmp2.get(k,0));
         }
         temp = cN.get(0,index);
         cN.put(0,index,cB.get(0,i));
         cB.put(0,i,temp);


         temp = x_N.get(index,0);
         x_N.put(index,0,x_B.get(i,0));
         x_B.put(i,0,temp);

        xB = B.getInverse()*b -B.getInverse()*N*xN;
         flag=true;
             for(int j=0;j<xN.getRows();j++)
             {
                 if(xN.get(j,0)<0)
                 {
                    flag=false;
                 }
             }
             for(int j=0;j<xB.getRows();j++)
             {
                 if(xB.get(j,0)<0)
                 {
                     flag=false;
                 }
             }
          if(flag)
          {
              return true;
          }
    }
    return false;
}

bool Simplex(Matrix<double> &B,Matrix<double> &N,Matrix<double> &b,Matrix<double> &xN,Matrix<double> &cN,Matrix<double> &cB,int iter,Matrix<double> & x_B,Matrix<double> &x_N)
{
    Matrix<double> xB = B.getInverse()*b -B.getInverse()*N*xN;
    Matrix<double> tmp = cN - cB*B.getInverse()*N;


    bool stop=true;

    for(int i =0;i<tmp.getColumns();i++)          //Stop criteria
    {
        if(tmp.get(0,i)>=0)
        {
            stop=false;
            break;
        }
    }
    iter=iter-1;
    if(stop||iter<=0)
    {
        for(int i=0;i<xN.getRows();i++)
        {
            if(xN.get(i,0)<0)
            {
                return false;
            }
        }
        for(int i=0;i<xB.getRows();i++)
        {
            if(xB.get(i,0)<0)
            {
                return false;
            }
        }
        return true;
    }

    int enter=-1;
    int leave=-1;
    double max=-999999;
    double min= 999999;
    for(int i =0;i<tmp.getColumns();i++)//Variable entering basis
    {
        if(tmp.get(0,i)>max)
        {
            max=tmp.get(0,i);
            enter=i;
        }
    }

    Matrix<double> tmp3 = B.getInverse()*N.getColumn(enter);
    Matrix<double> tmp2(xB.getRows(),1);
    for(int i=0;i<xB.getRows();i++)
    {
        tmp2.put(i,0,xB.get(i,0)/tmp3.get(i,0));
    }


    for(int i=0;i<tmp2.getRows();i++)////Variable leaving the basis
    {
        if(tmp2.get(i,0)<min&&tmp2.get(i,0)>=0)
        {
            min=tmp2.get(i,0);
            leave=i;
        }
    }

    tmp2 = B.getColumn(leave);//Swap the columns
    for(int i=0;i<B.getRows();i++)
    {
        B.put(i,leave,N.get(i,enter));
        N.put(i,enter,tmp2.get(i,0));
    }
    double temp = cN.get(0,enter);
    cN.put(0,enter,cB.get(0,leave));
    cB.put(0,leave,temp);


    temp = x_N.get(enter,0);
    x_N.put(enter,0,x_B.get(leave,0));
    x_B.put(leave,0,temp);

    return Simplex(B,N,b,xN,cN,cB,iter,x_B,x_N);
}



void print_matrix(const Matrix<double> & M)
{
  for (unsigned r = 0; r < M.getRows(); ++r) {
    for (unsigned c = 0; c < M.getColumns(); ++c) {
      cout << M.get(r, c) << " ";
    }
    cout << endl;
  }
}
