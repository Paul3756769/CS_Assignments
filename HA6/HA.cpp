#include "XY.cpp"
#include "MonteCarloLattice_new.cpp"

void XY_run(std::vector<double> L_Beta_U, int running_id = -1);
void debug(std::vector<double> L_Beta_U);	
void XY_runall(std::vector<std::vector<double>> L_Beta_Upd, int id,std::string filename);
void PottsW_run(std::string filename, int L, int dim, double beta, int q, int updates, bool init_hot = false, int I = 1, int measurement_id = 0);
void PottsW_runall(std::string filename, std::vector<std::vector<double>> L_dim_Q_Beta_updates, int measurement_id, double I = 1);
void HA11(int measurement_id,std::string filename = "HA11");
void HA12(int measurement_id,std::string filename = "HA12");

int main()
{
   
    //HA11(999, "data/HA11_test");
    HA12(995, "data/HA12");
    return 0;
}

void HA11(int measurement_id,std::string filename)
{
    std::vector<int> Ls = {8,24,60};
    int updates = std::pow(2,14)+5000;
    std::vector<std::vector<double>> L_Beta_Upd;


    for (int L:Ls)
    {
        for (double beta = 0; beta < 2.01; beta+=0.05)
        {
            L_Beta_Upd.push_back({(double)L,beta,(double)updates});
        }
        
    }
    XY_runall(L_Beta_Upd,measurement_id,filename);
    
}
  
void HA12(int measurement_id,std::string filename )
{
    std::vector<int> Ls = {8,16,32,64,128};
    int dim = 2;
    int q = 3;
    double I = 1;

    double beta = std::log(1+std::sqrt(q));
    int updates = std::pow(2,17)+5000;
    std::vector<std::vector<double>> L_dim_Q_Beta_updates;
    std::string Filename = filename+"_"+std::to_string(measurement_id);
    for (int L:Ls)
    {
        L_dim_Q_Beta_updates.push_back({(double)L,(double)dim,(double)q,beta,(double)updates});
    }
    PottsW_runall(Filename,L_dim_Q_Beta_updates,measurement_id,I);

}


void XY_run(std::vector<double> L_Beta_U, int running_id,std::string filename)
{
    double I = 1;
    double beta = L_Beta_U[1];
    int L = (int)L_Beta_U[0];
    int updates = (int)L_Beta_U[2];
    int measurement_id = (int)L_Beta_U[3];
    std::string running_id_str = "_"+std::to_string(running_id);
    if(running_id < 0)
        running_id_str = "";
    std::string Filename = filename+"_"+std::to_string(measurement_id)+running_id_str;

    XY_2D(Filename,L,beta,updates,I,measurement_id);
    std::cout << "Finished: " << Filename << std::endl;

}

void XY_runall(std::vector<std::vector<double>> L_Beta_Upd,int measurement_id,std::string filename)
{
    int running_id = 0;
    for (int i = 0; i < L_Beta_Upd.size(); i++)
    {
        L_Beta_Upd[i].push_back(measurement_id);
        XY_run(L_Beta_Upd[i], running_id, filename);
        running_id++;
    }
    
}

void PottsW_run(std::string filename, int L, int dim, double beta, int q, int updates, bool init_hot, int I, int measurement_id)
{
    std::set<int> images = {0};
    PottsWolffExperiment experiment(L,dim,q,beta,updates,images,init_hot,I);
    experiment.run_to_file_no_img(filename,measurement_id);
    std::cout << "Finished: " << filename << std::endl;
}
void PottsW_runall(std::string filename, std::vector<std::vector<double>> L_dim_Q_Beta_updates, int measurement_id, double I)
{
    int running_id = 0;
    for (int i = 0; i < L_dim_Q_Beta_updates.size(); i++)
    {
        std::vector<double> values = L_dim_Q_Beta_updates[i];
        int L = (int)values[0];
        int dim = (int)values[1];
        int q = (int)values[2];
        double beta = values[3];
        int updates = (int)values[4];
        std::string running_filename = filename+"_"+std::to_string(running_id);
        PottsW_run(running_filename,L,dim,beta,q,updates,true,I,measurement_id);
        running_id++;
    }
    
}