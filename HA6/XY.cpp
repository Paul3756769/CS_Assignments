#include <random>
#include <vector>
#include <map>	
#include <stack>
#include <fstream>
#include <iostream>
//#include <math.h>

const double PI = 3.14159265;//358979323846;
const int EXPERIMENT_ID = 5;

std::random_device rd;

std::mt19937 RND(rd());

void addToCluster(int index, double ising_angle,std::stack<int>* stack);
double flip(double angle,double ising_angle);
std::vector<double> init_lattice(int N);
double random_angle();
double random_index(int N);
double calc_energy(std::vector<double> spins,int L, int dim,double I = 1);
double calc_dE(std::vector<bool> spins, int index, int L, int dim);
double calc_magnetization(std::vector<double> spins);
std::stack<int> get_neighbors(int index,int L,int dimension);
int shift_one(int index,int dimension, int L, bool fwd = true);
bool aligned(double a,double ref);
void flip_spin(std::vector<double> *spins, int index, double ising_angle);
void XY_2D(std::string filename,int Length,double Beta, int Updates,double I=1, int id_value = 0);
bool criterion(double angle, double neighbor_angle, double ising_angle, double beta, double I=1);
double calc_p(double angle, double neighbor_angle, double ising_angle, double beta, double I);
void write_data_to_file(std::string filename, std::map<std::string, double> parameters, std::map<std::string,std::vector<double>> measurements, int id = 0);
double calc_m2(std::vector<double> spins);


void XY_2D(std::string Filename,int Length,double Beta, int Updates,double I, int id_value)
{
    const int dim = 2;
    const int L = Length;
    const int N = L*L;
    const double beta = Beta;
    int updates = Updates;
    std::string filename = Filename;
    
    std::map<std::string,std::vector<double>> timeseries;
    
    timeseries["E"] = std::vector<double>(updates,0);
    timeseries["M"] = std::vector<double>(updates,0);
    timeseries["m2"]= std::vector<double>(updates,0);
    timeseries["C"] = std::vector<double>(updates,0);

    std::vector<double> spins = init_lattice(N);

    double E = calc_energy(spins,L,dim,I);
    double M = calc_magnetization(spins);

    // build clusters
    for(int u = 0; u < updates; u++)
    {
        int cluster_size = 1;

        int start_index = random_index(N);
        double start_angle = spins[start_index];

        double ising_angle = random_angle();
        bool cluster_orientation = aligned(start_angle,ising_angle);


        std::stack<int> unchecked;
        unchecked.push(start_index);
        flip_spin(&spins,start_index,ising_angle);


        while(!unchecked.empty())
        {
            int index = unchecked.top();
            double angle = spins[index];
            unchecked.pop();

            for(int d=0;d<dim;d++)
            {
                for(bool fwd:{true,false})
                {
                    int neighbor = shift_one(index,d,L,fwd);
                    double neighbor_angle = spins[neighbor];

                    if (aligned(neighbor_angle,ising_angle) == cluster_orientation)
                    {
                        if(criterion(angle,neighbor_angle,ising_angle,beta,I))
                        {
                            unchecked.push(neighbor);
                            flip_spin(&spins,neighbor,ising_angle);
                            cluster_size++;
                        }
                    }

                }
            }
        }

        double new_E = calc_energy(spins,L,dim,I);
        double new_M = calc_magnetization(spins);
        double new_m2 = calc_m2(spins); 
        timeseries["E"][u] = new_E;
        timeseries["M"][u] = new_M;
        timeseries["C"][u] = cluster_size;
        timeseries["m2"][u] = new_m2;

    }


    // TODO File printing

        // TODO params
        std::map<std::string,double> params;
        params["L"] = L;
        params["beta"] = beta;
        params["I"] = I;
        params["number_of_updates"] = updates;
        params["dimension"] = dim;
        params["N"] = N;

        write_data_to_file(filename,params,timeseries,id_value);
}

bool criterion(double angle, double neighbor_angle, double ising_angle, double beta, double I)
{
    double p = calc_p(angle,neighbor_angle,ising_angle,beta,I);

    std::uniform_real_distribution<> dis(0, 1);
    return dis(RND) < p;
}
double calc_p(double angle, double neighbor_angle, double ising_angle, double beta, double I)
{
    double p = 1-exp(2.0*beta*I*cos(angle-ising_angle)*cos(neighbor_angle-ising_angle));
    return p;
}

void flip_spin(std::vector<double> *spins, int index, double ising_angle) {

  double *old_angle = &spins->at(index);

  double new_angle = flip(*old_angle, ising_angle);

  *old_angle = new_angle;
}

void addToCluster(int index, double ising_angle,std::stack<int>* stack)
{
        // TODO Delete here and in header
}

double correct_angle(double angle)
{
    double new_angle = angle;
    while (new_angle >= 2*PI || new_angle < 0)
    {
        if (new_angle >= 2*PI)
        {
            new_angle -= 2*PI;
        }
        else
        {
            new_angle += 2*PI;
        }
    }
    return new_angle; 
}

void correct_angle_new(double *angle)
{
    while (*angle >= 2 * PI || *angle < 0)
    {
        if (*angle >= 2 * PI)
        {
            *angle -= 2 * PI;
        }
        else
        {
            *angle += 2 * PI;
        }
    }
}

double flip(double angle,double ising_angle)
{
    double new_angle = 2*ising_angle-angle+PI;
    return correct_angle(new_angle);
}

std::vector<double> init_lattice(int N)
{
    std::vector<double> spins(N,0);
    for(int i=0;i<N;i++)
    {
        spins[i] = random_angle();
    }
    return spins;
}

double random_angle()
{
    std::uniform_real_distribution<> dis(0, 2*PI);
    return dis(RND);
}
double random_index(int N)
{
    std::uniform_int_distribution<> dis(0, N-1);
    return dis(RND);
}

double calc_energy(std::vector<double> spins,int L, int dim,double I)
{
    double E = 0;

    for(int i=0;i<spins.size();i++)
    {
        for(int d=0; d<dim;d++)
        {
            int neighbor_fwd = shift_one(i,d,L,true);
            E -= cos(spins[i]-spins[neighbor_fwd]);
        }
    }
    E*=I;

    return E;
}
double calc_dE(std::vector<bool> spins, int index, int L, int dim,double I) // TODO Funktion funktioniert nicht
{
    double dE = 0;
    for(int d=0;d<dim;d++)
    {
        for(bool fwd:{true,false})
        {
            int neighbor = shift_one(index,d,L,fwd);
            
        }
    }
    return 0;
}

double calc_magnetization(std::vector<double> spins) // TODO Definition von U. Wolff implementieren
{
    double M = 0;
    for(double spin:spins)
    {
        M += cos(spin);
    }
    return M;
}
double calc_m2(std::vector<double> spins) // TODO Definition von U. Wolff implementieren
{
    double M_x = 0;
    double M_y = 0;
    int V = spins.size();
    for(double spin:spins)
    {
        M_x += cos(spin);
        M_y += sin(spin); //TODO: Effizienz
    }
    double M2 = M_x*M_x+M_y*M_y;
    return M2/V/V;
}

std::stack<int> get_neighbors(int index,int L,int dimension)
{
    std::stack<int> neighbors;
    for(int d=0;d<dimension;d++)
    {
        neighbors.push(shift_one(index,d,L,true));
        neighbors.push(shift_one(index,d,L,false));
    } 
    return neighbors;
}

int shift_one(int index,int dimension, int L, bool fwd)
{
    int new_index = index;

        int test_index = index;

        // Fall: Index liegt am Rand
        test_index =  (int)((test_index % (int)pow(L,dimension+1))/pow(L,dimension));

        if (fwd)
        {
            new_index += pow(L,dimension);
            
            if (test_index==L-1)
            {
                new_index -= pow(L,dimension+1);
            }
            
        }
        else
        {
            new_index -= pow(L,dimension);
            
            if (test_index==0)
            {
                new_index += pow(L,dimension+1);
            }
            
        }
        return new_index;
}

void write_data_to_file(std::string filename, std::map<std::string, double> parameters,
                        std::map<std::string,std::vector<double>> measurements, int id)

{

  // Open file for writing
  std::ofstream outfile(filename+".dat");

  // Write header
  outfile << "Simulation (exp. id, measurement  id);"<<std::to_string(EXPERIMENT_ID)<<";" << std::to_string(id) << std::endl;

  // Write parameters
  bool first0 = true;
  for (auto& param : parameters) {
    if(!first0) outfile << ";"; else  first0 = false;
    outfile << param.first;
  }
  outfile << std::endl;
  first0 = true;
  for (auto& param : parameters) 
  {
    if(!first0) outfile << ";"; else first0 = false;
    outfile << param.second;
  }
  outfile << std::endl << std::endl;

  // Write image placeholders (the reading algorithm wants to see lattice images here)
 
    outfile<<0<<std::endl<<0<<std::endl<<std::endl;
  

  // Write column names
    first0 = true;
    std::vector<std::string> column_names;
    for (auto& quant : measurements) 
    {
        if(!first0)
        {
            outfile << ";"; 
        }
        else
        { 
            first0 = false;
        }
        outfile << quant.first;
        column_names.push_back(quant.first);
    }
    outfile << std::endl;
  // Write measurements
    int number_of_rows = measurements[column_names[0]].size(); 
    for (int row = 0; row < number_of_rows; row++)
    {
        bool first1 = true;
        for (auto& quant:measurements)
        {
            if(!first1) outfile << ";"; else  first1 = false;
            outfile << quant.second[row];
        }
        outfile << std::endl;
    }
    
  // Close file
  outfile.close();
}

bool aligned(double a,double ref)
{
    double A =correct_angle(a);
    double diff1 = std::abs(A-ref);
    bool dir = diff1<=PI/2 || diff1>3*PI/2;
    return dir;
}