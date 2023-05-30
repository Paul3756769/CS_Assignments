#include "MonteCarloLattice_new.cpp"



#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>


void HW3_P5(std::string filename);
void HW3_P6(std::string filename);
void HW4_P7(std::string filename);
void HW4_P8(std::string filename);
void HA5(std::string filename,int id);



int main()
{
    //HW3_P5("HW3_P5");
    //HW3_P6("HW3_P6");
    //HW3_P6("HW3_P6");
    //HW4_P7("HW4_P7");
    //HW4_P8("HW4_P8");

    int id = 9;
    HA5("data_HA5/HA5_A10_2_"+std::to_string(id),id);

    return 0;
}

void HW3_P5(std::string filename)
{
    int L[] = {4,8,16,32,64,128};
    //int L[] = {4};
    int dim = 2;
    double I = 1;
    bool hot = true;
    int sweeps = 10000;
    std::set<int> image_sweeps = {0};

    double critical_temp = 2/(log(1+sqrt(2)));

    //double beta[] = {1/critical_temp/(0.9),1/critical_temp/0.95,1/critical_temp,1/critical_temp/1.05,1/critical_temp/1.1};
    // 10 values
    double beta[] = {1/critical_temp/(0.9),1/critical_temp/(0.95),1/critical_temp/(0.975),1/critical_temp/(0.99),1/critical_temp/(0.995),1/critical_temp,1/critical_temp/(1.005),1/critical_temp/(1.01),1/critical_temp/(1.025),1/critical_temp/(1.05)};
    //double beta[] = {1/critical_temp};

    for(int l:L)
    {    
        int i = 1;
        for(double b:beta)
        {
            IsingWolffExperiment HW4(b,l,dim,I,hot,sweeps,image_sweeps);
            HW4.run_to_file(filename+"_"+std::to_string(l)+"_"+std::to_string(i));
            i++;
        }
    }
}
void HW3_P6(std::string filename)
{
    int L[] = {4,8,16};
    int dim = 2;
    double I = 1;
    bool hot = true;
    int sweeps = 10000;
    int q = 3;
    int beta = log(1+sqrt(q));
    std::set<int> image_sweeps = {0,100,1000,9999};

    for(int l:L)
    {    
        PottsWolffExperiment HW4(l,dim,q,beta,sweeps,image_sweeps,hot,I);
        HW4.run_to_file(filename+"_"+std::to_string(l));
    }
    

}

void HW4_P7(std::string filename)
{
    int L[] = {8,16,32,64};
    int dim = 2;
    double I = 1;
    bool hot = true;
    int sweeps = 100000;
    double beta_C = log(1+sqrt(2))/2;
    double Around[] = {0.95,1,1.05};
    std::set<int> image_sweeps = {0,100,1000,9999};
    for(int l:L)
    { 
        for(double around:Around)	
        {   

        double beta = beta_C*around;
        IsingWolffExperiment HW4(beta,l,dim,I,hot,sweeps,image_sweeps);
        int p = (int)(100*around);
        HW4.run_to_file(filename+"_"+std::to_string(l)+"_"+std::to_string(p));

        IsingMetropolisExperiment HW4_M(beta,l,dim,I,hot,sweeps,image_sweeps);
        HW4_M.run_to_file(filename+"_M_"+std::to_string(l)+"_"+std::to_string(p));
         
        }
    }
}

void HW4_P8(std::string filename)
{
    int L[] = {4,16,8};
    int dim = 2;
    double I = 1;
    bool hot = true;
    int sweeps = 100000;
    double beta_C = log(1+sqrt(2))/2;
    double Around[]  = {0.7,0.5,0.3,0.1,0.05,0.01};
    std::set<int> image_sweeps = {0,100,1000,9999};

    for(int l:L)
    {
      for (double around:Around)
      {
        double beta = beta_C*around;

        IsingWolffExperiment HW4(beta,l,dim,I,hot,sweeps,image_sweeps);
        HW4.run_to_file(filename+"_"+std::to_string(l)+"_"+std::to_string(beta));
    }
    }

}	
void HA5(std::string filename,int id)
{
    int L = 16;
    int dim = 2;
    double I = 1;
    bool hot = true;
    int sweeps = 5000+65536;
    double beta = log(1+sqrt(2))/2;

    IsingMetropolisExperiment HW5_M(beta,L,dim,I,hot,sweeps);
    HW5_M.run_to_file_no_img(filename,id);
}