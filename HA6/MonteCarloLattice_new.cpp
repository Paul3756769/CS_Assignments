#include <iostream>
#include <fstream>
#include <stdint.h>
#include <random>
#include <chrono>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <fstream>
#include <set>
#include <stack>


class Lattice
{
public:
    Lattice(int length, int dimension)
    {
        L = length;
        dim = dimension;
        N = pow(length,dimension);
        init_rnd();
    }
private:
    int L;
    int dim;
    int N;
      
    void init_rnd()
    {
        int sd = std::chrono::system_clock::now().time_since_epoch().count()%10000;
        rnd.seed(sd);
    }

protected:
    int shift_one(int index, int dimension, bool fwd = true)
    // Verschiebt den Index vorwärts/rückwärts in der angegebenen Richtung (dimension [0,1,2...])
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
    int select_random_index()
    {
        std::uniform_int_distribution<int> dist(0,N-1);
        return dist(rnd);
    }
    
    int get_L()
    {
        return L;
    }
    int get_dim()
    {
        return dim;
    }
    int get_N()
    {
        return N;
    }

     std::mt19937 rnd;

    std::vector<int> get_neighbors(int index)
    {
        std::vector<int> neighbors;
        for(int i=0;i<dim;i++)
        {
            neighbors.push_back(shift_one(index,i,true));
            neighbors.push_back(shift_one(index,i,false));
        }
        return neighbors;
    }
};
class PottsModel : public Lattice
{
    public:
    PottsModel(int length, int dimension, int q, bool init_hot = false, double coupling_constant = 1) : Lattice(length,dimension)
    {
        init_model(length,dimension,q,init_hot,coupling_constant);
    }

    protected:
    int select_random_value()
    {
        std::uniform_int_distribution<int> dist(0,get_q()-1);
        return dist(rnd);
    }
    int select_random_value(int q0)
    {
        std::uniform_int_distribution<int> dist(0,get_q()-1);
        int q_new = dist(rnd);
        while (q_new == q0)
        {
            std::uniform_int_distribution<int> dist2(0,get_q()-1);
            q_new = dist2(rnd);
        }
        
        if (q_new == q0)
        {
            throw std::runtime_error("Error: q_new == q0");
        }
        
        return q_new;
    }

    double get_E()
    {
        return E;
    }
    int get_t()
    {
        return time;
    }
    protected:

    void set_t(int t=0)
    {
        time = t;
    }
    void inc_t()
    {
        time++;
    }

    void set_E(double energy)
    {
        E = energy;
    }
    void update_E()
    {
        set_E(calc_E());
    }
    double calc_E()
    {
        double energy = 0;
        for(int i=0;i<get_N();i++)
        {
            std::vector<int> neighbors = get_neighbors(i);
            double cc = get_I();
            int this_value = get_Value(i);
            for(int neighbor:neighbors)
            {
                if(this_value == get_Value(neighbor)) 
                    energy -= cc;
            }
        }
        return energy;
    }

    int get_q()
    {
        return q;
    }
    void set_q(int q_new)
    {
        q = q_new;
    }
    double get_I()
    {
        return I;
    }
    int get_Value(int index)
    {
        return lattice[index];
    }
    void set_Value(int index, int value)
    {
        if(verify_q(value))   lattice[index] = value;
    }
    bool verify_q(int q_new)
    {
        if(q_new<0 || q_new>get_q()) return false;
        else return true;
    }



    private:


    double E;
    int time;

    void init_model(int length, int dimension, int Q, bool init_hot = false, double coupling_constant = 1)
    {
        
        set_q(Q);
        I = coupling_constant;
        E = 0;
        time = 0;
        init_lattice(init_hot);
    }
    int q;
    double I;
    std::vector<int> lattice;
    void init_lattice(bool init_hot = false)
    {
        lattice.resize(get_N());
        if(init_hot)
        {
            std::uniform_int_distribution<int> dist(0,get_q()-1);
            for(int i=0;i<get_N();i++)
            {
                //lattice[i] = dist(rnd);
                set_Value(i,dist(rnd));
            }
        }
        else
        {
            for(int i=0;i<get_N();i++)
            {
                //lattice[i] = 0;
                set_Value(i,0);
            }
        }
        
    }
};
class PottsWolff : public PottsModel
{
    public:
    PottsWolff(int length, int dimension, double beta, int q, bool init_hot = false, double coupling_constant = 1) : PottsModel(length,dimension,q,init_hot,coupling_constant)
    {
        set_beta(beta);
    }
    protected:
    void set_beta(double beta)
    {
        Beta = beta;
        define_p();
    }
    double get_beta()
    {
        return Beta;
    }
     void WolffSweep(bool justOne = false)
    {
        

        
        
        if (justOne)
        {
            int actual_cluster_size = WolffUpdate();
            total_cluster_size += actual_cluster_size;
            total_number_of_clusters++;
            current_cluster_size = (double)actual_cluster_size;

        }
        else
        {
            int accumulated_cluster_size = 0;
            int N = get_N();
            int actual_cluster_size = 0;
            int numberofclustersinsweep = 0;
            while (accumulated_cluster_size + 0.5*actual_cluster_size < N)
            {
                actual_cluster_size = WolffUpdate();
                accumulated_cluster_size += actual_cluster_size;
                total_number_of_clusters++;
                numberofclustersinsweep++;
            }
            total_cluster_size += accumulated_cluster_size;
            current_cluster_size = (double)accumulated_cluster_size/(double)numberofclustersinsweep;
        }
        
        
        update_E();
        inc_t();
    }
    protected:
    double get_current_cluster_size()
    {
        return current_cluster_size;
    }
    double get_average_cluster_size()
    {
        return (double)total_cluster_size/(double)total_number_of_clusters;
    }

    private:

    int total_cluster_size = 0;
    int total_number_of_clusters = 0;
    int current_cluster_size = 0;


    double Beta;
    
    double p;
    void define_p()
    {
        p = 1 - exp(-get_beta()*get_I());
    }

    int WolffUpdate()
    {
        std::stack<int> indices_toBeChecked;
        int starting_index = select_random_index();
        int q0 = get_Value(starting_index);
        int q_new = select_random_value();
        while (q0 == q_new)
        {
            q_new = select_random_value();
        }
        

        int cluster_size = 1;

        set_Value(starting_index,q_new);
        indices_toBeChecked.push(starting_index);
        while (!indices_toBeChecked.empty())
        {
            int current_index = indices_toBeChecked.top();
            indices_toBeChecked.pop();
            std::vector<int> neighbors = get_neighbors(current_index);
            for (int neighbor:neighbors)
            {
                    if (Criterion(neighbor,q0))
                    {
                        set_Value(neighbor,q_new);
                        indices_toBeChecked.push(neighbor);
                        cluster_size++;
                    }
             
            }
            // debug
            //std::cout << "cluster size: " << cluster_size <<";"<< indices_toBeChecked.size()<< std::endl;
            
        }
        int test = cluster_size;

        return cluster_size;
    }

   
    

    bool Criterion(int index, int q0)
    {
        if (get_Value(index) == q0)
        {
            std::uniform_real_distribution<double> dist(0,1);
            if (dist(rnd) < p) return true;
            else return false;
        }
        else return false;
    }

    protected:
    int select_random_index()
    {
        std::uniform_int_distribution<int> dist(0,get_N()-1);
        return dist(rnd);
    }

};
class Ising
{
public:
    Ising(int length, int dimension,double coupling_constant=1,bool init_temp=false, double spinvalue=1)
    {
        init_model(length,dimension,coupling_constant,init_temp,spinvalue);
    }
   
private:
    //Gitter
    std::vector<int> lattice;
    void init_lattice(bool init_temp=false)
    {
        lattice.resize(N);
        if(init_temp)
        {
            std::uniform_int_distribution<int> dist(0,1);
            for(int i=0;i<N;i++)
            {
                //lattice[i] = dist(rnd);
                int value = 2*dist(rnd)-1;
                set_Spin(i,value);
            }
        }
        else
        {
            for(int i=0;i<N;i++)
            {
                //lattice[i] = 0;
                set_Spin(i,0);
            }
        }
    }
    // Parameter
    //double beta; //inverse Temperatur
    int L;      //Seitenlänge des Gitters
    int dim;    //Dimension des Gitters (1,2,3...)
    int N;      //Anzahl der Gitterelemente
    double I=1;  //Kopplungskonstante
    double abs_spin = 1; //Betrag des Spins

    //Physikalische Größen
    double E; //Energie
    double M; //Magnetisierung
    int t;      //Zeit (Sweeps)

protected:
    // Hilfsfunktionen
    void init_model(int length, int dimension,double coupling_constant=1,bool init_temp=false, double spinvalue=1)
    {
        // Parameterwerte setzen
        L = length;
        dim = dimension;
        I = coupling_constant;
        N = pow(L,dim);
        // Gitter initialisieren
        init_lattice(init_temp);
        // Physikalische Größen berechnen
        E = calc_E();
        M = calc_M();
        t = 0;
        abs_spin = spinvalue;
    }
    int shift_one(int index, int dimension, bool fwd = true)
    // Verschiebt den Index vorwärts/rückwärts in der angegebenen Richtung (dimension [0,1,2...])
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
    void flip(int index)
    {
        //lattice[index] = !lattice[index];
        if (get_Spin(index)==1) set_Spin(index,-1);
        else set_Spin(index,1);
    }	
    double calc_E()
    // Energie aus aktuellem Gitter berechnen
    {
        double E = 0;
        for(int i=0;i<N;i++)
        {
            double spin = get_Spin(i);
            for(int d=0;d<dim;d++)
            {
                int fw = shift_one(i,d);
                E += spin*get_Spin(fw);
            }
        }

        return E*I*(-1);
    }
    double calc_M()
    // Magnetisierung aus aktuellem Gitter berechnen
    {
        double M = 0;
        for(int i=0;i<N;i++)
        {
            M += get_Spin(i);
        }
        return M;
    }
    double get_dE(int index)
    {
        double neighborspins = 0;
        double spin = get_Spin(index);
        for(int d=0;d<dim;d++)
        {
            int fw = shift_one(index,d,true);
            int bw = shift_one(index,d,false);
            neighborspins += get_Spin(fw) + get_Spin(bw);
            
        }
        return  2*I*spin*neighborspins;
    }

    //Zufallszahlengenerator
    
    std::mt19937 rnd // Mersenne twister, Seed = 0 oder Systemzeit
    {
        0
        //std::chrono::high_resolution_clock::now().time_since_epoch().count()  % 10000
    };
    int select_random_index()
    {
        std::uniform_int_distribution<int> dist(0,N-1);
        return dist(rnd);
    }

    
    
    std::vector<int> get_lattice()
    {
        return lattice;
    }
    //Abrufen der Parameter
    double get_abs_spin()
    {
        return abs_spin;
    } 
    int get_L()
    {
        return L;
    }
    int get_dim()
    {
        return dim;
    }
    int get_N()
    {
        return N;
    }
    double get_I()
    {
        return I;
    }
    double get_E()
    {
        return E;
    }
    void update_E()
    {
        E = calc_E();
    }
    void set_E(double new_E)
    {
        E = new_E;
    }
    void add_dE(double dE)
    {
        E += dE;
    }
    double get_M()
    {
        return M;
    }
    void update_M()
    {
        M  = calc_M();
    }
    void set_M(double new_M)
    {
        M = new_M;
    }
    void add_dM(double dM)
    {
        M += dM;
    }
    int get_t()
    {
        return t;
    }
    void add_t(int dt=1)
    {
        t += dt;
    }
    void reset_t(int new_t=0)
    {
        t = new_t;
    }
    double get_Spin(int index)
    {
        if(lattice[index]==1) return abs_spin;
        else return -abs_spin;
    }
    void set_Spin(int index, int new_spin)
    {
        if(new_spin == 1) lattice[index] = 1;
        else if(new_spin ==-1) lattice[index] = -1;
        //lattice[index] = new_spin;
    }
};
class IsingMetropolis : public Ising
{
public:
    //Konstruktor
    IsingMetropolis(double beta, int length,int dimension,double coupling_constant=1,bool init_temp=false):Ising(length,dimension,coupling_constant,init_temp)
    {
        set_Beta(beta);
        exp_table = init_exp_table();
    }
    void MetropolisSweep()
    {
        int N_ = get_N();

        for(int i=0;i<N_;i++)
        {
            Metropolis_Update();
        }
        
        add_t();
    }

protected:
    //Physikalische Größen
    void set_Beta(double new_Beta)
    {
        Beta = new_Beta;
    }
    double get_Beta()
    {
        return Beta;
    }

private:
//private
    double Beta; //inverse Temperatur
    //Hilfsfunktionen
    std::map<double,double> exp_table; // Wertetabelle füt die Gibbs-Verteilung
    std::map<double,double> init_exp_table()
    {
        double Coupling_C = get_I();
        int Dimensions = get_dim();
        double s = get_abs_spin();
        //double s2 = s*s;                    //Wert eines Spins
        std::map<double,double> hashmap;
        //hashmap[0] = 1;                     // exp 0 = 1
        for(int dddd=0;dddd<=get_dim();dddd++)      //Dimensions+1 ist die Anzahl der möglichen Werte für dE
        {
            double ds = 2*dddd*s*s;
           // double E_current = 2*d*s2;
            double dE = 2*Coupling_C*ds;
            //std::cout << "dE = " << dE << std::endl;
            hashmap[dE] = exp(-dE*get_Beta());
            hashmap[-dE] = exp(dE*get_Beta());    // eigentlich überflüssig, da >1
        }
        return hashmap;     
    }
    double GibbsDistr(double dE)
    //Raum für Fehlerbehandlung
    {
        double p = exp_table[dE]; 
        return p;
    }
    bool Metropolis_Criterion(double dE)
    {
        if (dE<0)
        {
            return true;
        }
        else
        {
            std::uniform_real_distribution<double> dist(0,1);
            double r = dist(rnd);

            if (r<GibbsDistr(dE))   return true;
        
            else                    return false;
        }
    }
    void Metropolis_Update()
    {
        int index = select_random_index();
        double dE = get_dE(index);
        if (Metropolis_Criterion(dE))
        {
            flip(index);
            

            add_dE(dE);
            
            double dM = 2*get_Spin(index);
            add_dM(dM);
        }
    }
};
class IsingWolff : public Ising
{
public: 
    IsingWolff(double beta, int length,int dimension,double coupling_constant=1,bool init_temp=false):Ising(length,dimension,coupling_constant,init_temp)
    {
        set_Beta(beta);
    }
    void WolffSweep(bool justOne = false)
    {
        Clear_Stack();

        
        
        if (justOne)
        {
            int actual_cluster_size = Wolff_Update();
            total_cluster_size += actual_cluster_size;
            total_number_of_clusters++;
            current_cluster_size = (double)actual_cluster_size;

        }
        else
        {
            int accumulated_cluster_size = 0;
            int N = get_N();
            int actual_cluster_size = 0;
            int numberofclustersinsweep = 0;
            while (accumulated_cluster_size + 0.5*actual_cluster_size < N)
            {
                actual_cluster_size = Wolff_Update();
                accumulated_cluster_size += actual_cluster_size;
                total_number_of_clusters++;
                numberofclustersinsweep++;
            }
            total_cluster_size += accumulated_cluster_size;
            current_cluster_size = (double)accumulated_cluster_size/(double)numberofclustersinsweep;
        }

        update_M();
        update_E();
        add_t();
    }
    
protected:
    //Physikalische Größen
    void set_Beta(double new_Beta)
    {
        Beta = new_Beta;
        define_binding_probability(new_Beta);
        define_SweepSize(new_Beta);
    }
    double get_Beta()
    {
        return Beta;
    }
    double get_total_average_cluster_size()
    {
        return (double)total_cluster_size/(double)total_number_of_clusters;
    }
    double get_current_cluster_size()
    {
        return current_cluster_size;
    }
private:   
    double Beta; //inverse Temperatur
    double binding_probability; //Wahrscheinlichkeit, dass ein Spin gebunden wird
    double SweepSize;
    
    int total_cluster_size= 0;
    int total_number_of_clusters= 0;
    

    double current_cluster_size = 0;
    


    std::stack<int> Stack;
    

    //Hilfsfunktionen
    double get_binding_probability()
    {
        return binding_probability;
    }
    void define_binding_probability(double beta)
    {
        binding_probability = 1-exp(-2*beta*get_I());
    }
    void define_SweepSize(double beta)
    {
        SweepSize = 1; //get_N();
    }
    double get_SweepSize()
    {
        return SweepSize;
    }
    int Wolff_Update() // Build Cluster
    {
        int starting_index = select_random_index();
        addToCluster(starting_index);
        int cluster_size = 1;
        while (StackNotEmpty())
        {
            int stacked_index = take_index_from_stack();
            for (int d = 0; d < get_dim(); d++)
            {
                for(bool b : {true,false})
                {
                    int neighbor_index = shift_one(stacked_index,d,b);
                    if (BindingCriterion(stacked_index,neighbor_index))
                        {
                            addToCluster(neighbor_index);
                            cluster_size++;
                        }
                }
            }
        }
        return cluster_size;
    }
    
    bool BindingCriterion(int index1, int index2)
    {
        if (get_Spin(index1)!=get_Spin(index2))
        {
            std::uniform_real_distribution<double> dist(0,1);
            double r = dist(rnd);
            if (r<get_binding_probability()) return true;
            else return false;
        }
        else return false;
    }
    void addToCluster(int index)
    {
        Stack.push(index);
        flip(index);
    }
    bool StackNotEmpty()
    {
        return !Stack.empty();
    }
    int take_index_from_stack()
    {
        int index = Stack.top();
        Stack.pop();
        return index;
    }
    void Clear_Stack()
    {
        while (StackNotEmpty())
        {
            Stack.pop();
        }
    }

    
};
class MonteCarloExperiment
{
        
public:
    MonteCarloExperiment(int number_of_sweeps_, std::set<int> sweepnumbers_of_lattice_images_)
    {
        number_of_sweeps = number_of_sweeps_;
        numbers_of_lattice_images = sweepnumbers_of_lattice_images_;
        E_timeseries.resize(number_of_sweeps+1);
        M_timeseries.resize(number_of_sweeps+1);
        t_timeseries.resize(number_of_sweeps+1);

        init_measurements_new(quantity_names);
        
        init_default_parameters();
        

        
    }
    
    void run_to_file(std::string Filename)
    {
        generate_data();
        std::vector<std::string> Quantity_names = quantity_names;
        std::vector<std::vector<double>> Measurements = measurements;
        std::map<std::string,double> params = parameters;
        std::map<std::string,std::vector<double>> Measurements_new = measurements_new;
        std::set<int> image_times = numbers_of_lattice_images;
        std::map<int,std::vector<int>> Lattice_images = lattice_images_new;

        //write_data_to_file
        write_data_to_file(Filename,params,Measurements_new,Lattice_images);
    }

    void run_to_file_no_img(std::string filename, int id)
    {
        generate_data();
        std::map<std::string,double> params = parameters;
        std::map<std::string,std::vector<double>> Measurements_new = measurements_new;
        
        //write_data_to_file
        write_data_to_file_no_img(filename,params,Measurements_new,id);
    }


int EXPERIMENT_ID = 6;
void write_data_to_file_no_img(std::string filename, std::map<std::string, double> parameters,
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
}

    void run(std::string Filename)
    {
        generate_data();
       
        std::ofstream file(Filename+".csv");	
        for (int i = 0; i < quantity_names.size(); i++)
        {
            file << quantity_names[i]<<";";
        }
        file << "\n";
        for (int i = 0; i < measurements[0].size(); i++)
        {
            for (int j = 0; j < quantity_names.size(); j++)
            {
                int k = measured_quantities[quantity_names[j]];
                file << measurements[k][i]<<";";
            }
            file << "\n";
        }
        file.close();
        std::ofstream file2(Filename+"_lattice_images.csv");
        for (int i = 0; i < lattice_images.size(); i++)
        {
            file2 << lattice_to_string(lattice_images[i])<<"\n";
        }
        
        
    }
    void run_to_console()
    {
        generate_data();
            for (int i = 0; i < quantity_names.size(); i++)
            {
                std::cout << quantity_names[i]<<";";
            }
            std::cout << "\n";
            for (int i = 0; i < measurements[0].size(); i++)
            {
                for (int j = 0; j < quantity_names.size(); j++)
                {
                    int k = measured_quantities[quantity_names[j]];
                    std::cout << measurements[k][i]<<";";
                }
                std::cout << "\n";
            }
        std::cout << "\n";
        print_lattice_images();
    }

protected:
    std::vector<double> E_timeseries;
    std::vector<double> M_timeseries;
    std::vector<int> t_timeseries;
    
    
    std::vector<std::vector<double>> measurements;
    std::map<std::string,std::vector<double>> measurements_new;
    std::map<std::string,int> measured_quantities; //maps names of quantities to their index in measurements
    std::vector<std::string> quantity_names;
    std::vector<std::vector<int>> lattice_images;
    std::set<int> numbers_of_lattice_images;
    std::map<int,std::vector<int>> lattice_images_new;
    std::map<std::string,double> parameters;
    void init_default_parameters() //initializes parameters
    {
        parameters["number_of_sweeps"] = number_of_sweeps;
    }



    virtual void init_Measurements() = 0; //initializes measurements and measured_quantities

    //one measurement will be performed after each sweep

    int number_of_sweeps; // = number of measurements -1

    

    void init_measurements_new(std::vector<std::string> Quantity_names)
    {
        for (int i = 0; i < Quantity_names.size(); i++)
        {
            measurements_new[Quantity_names[i]] = std::vector<double>(number_of_sweeps+1);
        }
        
    }
    
    std::vector<double> get_E_timeseries()
    {
        return E_timeseries;
    }
    std::vector<double> get_M_timeseries()
    {
        return M_timeseries;
    }
    std::vector<int> get_t_timeseries()
    {
        return t_timeseries;
    }
    
    void generate_data()
    {
        measure(0);
        for(int i=0;i<number_of_sweeps;i++)
        {
            sweep();
            measure(i+1);
        }
    }

    virtual void measure(int measurement_index) = 0;

    virtual void sweep() = 0;

    void print_lattice_images()
    {
        for (int i = 0; i < lattice_images.size(); i++)
        {
                std::cout << lattice_to_string(lattice_images[i]) << "\n";
        }
    }
    
     std::string lattice_to_string(std::vector<int> lattice)
    {
        std::string s = "";
        for (int i = 0; i < lattice.size(); i++)
        {
            s += std::to_string(lattice[i]);           
        }
        return s;
    }
    
/*void write_data_to_file(std::string filename, std::map<std::string, double> parameters,
                        std::vector<std::vector<double>> measurements,
                        std::map<int, std::vector<int>> images)

// Based on a draft created by ChatGPT                        
{

  // Open file for writing
  std::ofstream outfile(filename);

  // Write header
  outfile << "Simulation data file" << std::endl;

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

  // Write images
  if (!images.empty()) 
  {
     first0 = true;
    for (auto& image : images) {
      if(!first0) outfile << ";"; else first0 = false;
      outfile << image.first;
    }
    outfile << std::endl;


    first0 = true;
    for (auto& image : images) {
      if(!first0) outfile << ";"; else first0 = false;
      std::vector<int> image_sequence = image.second;
      bool first1 = true;
      for (auto& val : image_sequence) {
        if(!first1) outfile << "_"; else first1 = false;
        outfile << val;
      
      }
    }
    outfile << std::endl << std::endl;
  }

  // Write column names
  outfile << "t;E;M"<< std::endl;

  // Write measurements
  for (auto& row : measurements) 
  {
    bool first1 = true;
    for (auto& val : row) {
      if(!first1) outfile << ";"; else  first1 = false;
      outfile << val;

    }
    outfile << std::endl;
  }

  // Close file
  outfile.close();
}*/

void write_data_to_file(std::string filename, std::map<std::string, double> parameters,
                        std::map<std::string,std::vector<double>> measurements,
                        std::map<int, std::vector<int>> images)

// Based on a draft created by ChatGPT                        
{

  // Open file for writing
  std::ofstream outfile(filename+".dat");

  // Write header
  outfile << "Simulation data file" << std::endl;

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

  // Write images
  if (!images.empty()) 
  {
    first0 = true;
    for (auto& image : images) 
    {
      if(!first0) outfile << ";"; else first0 = false;
      outfile << image.first;
    }
    outfile << std::endl;


    first0 = true;
    for (auto& image : images) 
    {
      if(!first0) outfile << ";"; else first0 = false;
      std::vector<int> image_sequence = image.second;
      bool first1 = true;
      for (auto& val : image_sequence) {
        if(!first1) outfile << "_"; else first1 = false;
        outfile << val;
      
      }
    }
    outfile << std::endl << std::endl;
  }
  else
  {
        outfile<<0<<std::endl<<0<<std::endl<<std::endl;
  }

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

};
class IsingExperiment : public MonteCarloExperiment, public Ising
{
public: 
    IsingExperiment(int length,int dimension, int number_of_sweeps_, bool init_hot = true,std::set<int> sweepnumbers_of_lattice_images_={},double coupling_constant_=1, double spinvalue = 1):Ising(length,dimension,coupling_constant_,init_hot,spinvalue),MonteCarloExperiment(number_of_sweeps_,sweepnumbers_of_lattice_images_)
    {
        init_Measurements();
        std::vector<std::string> Quantity_names = {"E","M","t"};
        init_measurements_new(Quantity_names);

    }

private:
    void init_Measurements()
    {
        quantity_names = {"E","M","t"};

        measurements.resize(quantity_names.size());
        for(int i=0;i<quantity_names.size();i++)
        {
            measured_quantities[quantity_names[i]] = i;
            measurements[i].resize(number_of_sweeps+1);

        }
    }
    void measure(int measurement_index)
    {
        E_timeseries[measurement_index] = get_E();
        M_timeseries[measurement_index] = get_M();
        t_timeseries[measurement_index] = get_t();

        measurements[measured_quantities["E"]][measurement_index] = get_E();
        measurements[measured_quantities["M"]][measurement_index] = get_M();
        measurements[measured_quantities["t"]][measurement_index] = get_t();

        measurements_new["E"][measurement_index] = get_E();
        measurements_new["M"][measurement_index] = get_M();
        measurements_new["t"][measurement_index] = get_t();
    }
};
class IsingMetropolisExperiment : public MonteCarloExperiment, public IsingMetropolis
{
public:
    IsingMetropolisExperiment(double beta, int length,int dimension,double coupling_constant=1,bool init_temp=false,int number_of_sweeps_=1000, std::set<int> sweepnumbers_of_lattice_images_={0}):IsingMetropolis(beta,length,dimension,coupling_constant,init_temp),MonteCarloExperiment(number_of_sweeps_,sweepnumbers_of_lattice_images_)
    {
        
        parameters["beta"] = beta;
        parameters["length"] = length;
        parameters["dimension"] = dimension;
        parameters["coupling_constant"] = coupling_constant;
        parameters["init_temp"] = init_temp;
        parameters["abs_spin"] = get_abs_spin();
        init_Measurements();

        std::vector<std::string> Quantity_names = {"E","M","t"};
        init_measurements_new(Quantity_names);


    }
    void measure(int measurement_index)
    {
        E_timeseries[measurement_index] = get_E();
        M_timeseries[measurement_index] = get_M();
        t_timeseries[measurement_index] = get_t();

        measurements[measured_quantities["E"]][measurement_index] = get_E();
        measurements[measured_quantities["M"]][measurement_index] = get_M();
        measurements[measured_quantities["t"]][measurement_index] = get_t();

        if(numbers_of_lattice_images.find(measurement_index-1)!=numbers_of_lattice_images.end())
        {
            lattice_images.push_back(get_lattice());
        }
        if(numbers_of_lattice_images.find(measurement_index-1)!=numbers_of_lattice_images.end())
        {
            lattice_images_new[measurement_index-1] = get_lattice();
        }

        measurements_new["E"][measurement_index] = get_E();
        measurements_new["M"][measurement_index] = get_M();
        measurements_new["t"][measurement_index] = get_t();
      
        

    }
    void sweep() {MetropolisSweep();};
    void init_Measurements()
    {
        quantity_names = {"E","M","t"};

        measurements.resize(quantity_names.size());
        for(int i=0;i<quantity_names.size();i++)
        {
            measured_quantities[quantity_names[i]] = i;
            measurements[i].resize(number_of_sweeps+1);

        }

       
    }
};
class IsingWolffExperiment : public MonteCarloExperiment, public IsingWolff
{
    public:
    IsingWolffExperiment(double beta, int length,int dimension,double coupling_constant=1,bool init_temp=false,int number_of_sweeps_=1000, std::set<int> sweepnumbers_of_lattice_images_={0}):IsingWolff(beta,length,dimension,coupling_constant,init_temp),MonteCarloExperiment(number_of_sweeps_,sweepnumbers_of_lattice_images_)
    {
        init_Measurements();
        std::vector<std::string> Quantity_names = {"E","M","t","C"};
        init_measurements_new(Quantity_names);

         
        parameters["beta"] = beta;
        parameters["length"] = length;
        parameters["dimension"] = dimension;
        parameters["coupling_constant"] = coupling_constant;
        parameters["init_temp"] = init_temp;
        parameters["abs_spin"] = get_abs_spin();
    }
    void measure(int measurement_index)
    {
        measurements[measured_quantities["E"]][measurement_index] = get_E();
        measurements[measured_quantities["M"]][measurement_index] = get_M();
        measurements[measured_quantities["t"]][measurement_index] = get_t();


        if(numbers_of_lattice_images.find(measurement_index-1)!=numbers_of_lattice_images.end())
        {
            lattice_images.push_back(get_lattice());
        }

        if(numbers_of_lattice_images.find(measurement_index-1)!=numbers_of_lattice_images.end())
        {
            lattice_images_new[measurement_index-1] = get_lattice();
        }

        measurements_new["E"][measurement_index] = get_E();
        measurements_new["M"][measurement_index] = get_M();
        measurements_new["t"][measurement_index] = get_t();
        measurements_new["C"][measurement_index] = get_current_cluster_size();
      
        

    }
    void sweep() {WolffSweep(true);};
    void init_Measurements()
    {
        quantity_names = {"E","M","t"};

        measurements.resize(quantity_names.size());
        for(int i=0;i<quantity_names.size();i++)
        {
            measured_quantities[quantity_names[i]] = i;
            measurements[i].resize(number_of_sweeps+1);
        }

       
    }
};
class PottsWolffExperiment : public MonteCarloExperiment, public PottsWolff
{
    public:
    PottsWolffExperiment(int length, int dimension,int q, double beta, int number_of_sweeps_, std::set<int> sweepnumbers_of_lattice_images_, bool init_hot = false, double coupling_constant = 1) : PottsWolff(length,dimension,beta,q,init_hot,coupling_constant), MonteCarloExperiment(number_of_sweeps_,sweepnumbers_of_lattice_images_)
    {
        init_Measurements();
        std::vector<std::string> Quantity_names = {"E","t","C"};
        init_measurements_new(Quantity_names);

         
        parameters["beta"] = beta;
        parameters["length"] = length;
        parameters["dimension"] = dimension;
        parameters["coupling_constant"] = coupling_constant;
        parameters["init_temp"] = init_hot;
        parameters["q"] = q;
    }

   protected:
    void sweep()
    {
        WolffSweep(true);
    }

   private:



   void init_Measurements()
    {
        quantity_names = {"E","t","C"};

        measurements.resize(quantity_names.size());
        for(int i=0;i<quantity_names.size();i++)
        {
            measured_quantities[quantity_names[i]] = i;
            measurements[i].resize(number_of_sweeps+1);

        }
    }
    void measure(int measurement_index)
    {
        E_timeseries[measurement_index] = get_E();
       
        t_timeseries[measurement_index] = get_t();

        measurements[measured_quantities["E"]][measurement_index] = get_E();

        measurements[measured_quantities["t"]][measurement_index] = get_t();

       measurements_new["E"][measurement_index] = get_E();
       measurements_new["t"][measurement_index] = get_t();
        measurements_new["C"][measurement_index] = get_current_cluster_size();

    /*
    if(numbers_of_lattice_images.find(measurement_index-1)!=numbers_of_lattice_images.end())
        {
            lattice_images.push_back(get_lattice());
        }

        if(numbers_of_lattice_images.find(measurement_index-1)!=numbers_of_lattice_images.end())
        {
            lattice_images_new[measurement_index-1] = get_lattice();
        }
    */
    }
};

