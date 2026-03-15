#include <iostream>
#include <string>
#include <utility>  // For std::move
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>

int NUMBER_OF_PARTICLES = 1000; //param by default = ...
int TIME_STEPS_PER_DELTA = 365; //param by default = ...
double EPSILON = 1e-5; //0.00001 put inside calibrater-private: static constexpr double EPSILON = 1e-9;

double LA = 1;


double f0(double T) {
    return 0.12; // Market instantaneous forward rate
}

double P0(double T) {
    return std::exp(T * -0.12); // Market ZCB
}

double mm_part(double T, double K) {
    return 0.01; // here i have the correct market price
}

//so strikes can NOT be modified inside (ssince we pass reference!!) -> write const
double get_vol_from_smile(const std::vector<double>& strikes, const std::vector<double>& vols, double s) { // linear interpolation, wrt volatilities : Vol_value (K_grid)
        if (s <= strikes.front()) return vols.front();
        if (s >= strikes.back()) return vols.back();
        auto it = std::lower_bound(strikes.begin(), strikes.end(), s);
        size_t j = it - strikes.begin();
        size_t i = j - 1;
        double w = (s - strikes[i]) / (strikes[j] - strikes[i]);
        return vols[i] + w*(vols[j] - vols[i]);
} 


int main() {
    auto start = std::chrono::high_resolution_clock::now();
    auto G = [](double t1, double t2)
    {
        return (1 - std::exp(-LA*(t2 - t1))) / LA;
    };

 
    //initialization (allocate memory once))
    double xMax = 10;
    double nx = 10;
    double xMin = 0;

    double yMax;  //Get from CDF
    double ny = 30;
    double yMin;  //Get from CDF 

    double dx = (xMax - xMin) / nx; //OR (xMax - xMin) / (nx - 1)
    double dy;

    //use statically allocated vectors
    std::vector<std::vector<std::pair<double,double>>> grid(nx, std::vector<std::pair<double,double>>(ny));
    std::vector<std::vector<double>> vol_surface(nx, std::vector<double>(ny));
    std::vector<double> prev_Strike_grid(ny);
    std::vector<double> prev_Vol_values(ny);

    double T;  // build sigma(T, K) surface
    double K; 
    double term;
    double particle_dt = (xMax - xMin) / (nx * TIME_STEPS_PER_DELTA);
    double prev_s;
    double prev_dsdx ;
    int NumberFromFixLegSchedule;
    double ti;


    //---------------- INITIALIZE 
    term = xMax - EPSILON; //for co-term meturity or term = xMax
    //define yMax, yMin
    yMin = 0.1; //Get from CDF using 
    yMax = 1; //Get from CDF
    dy = (yMax - yMin) / ny;  //not including right limit

    for (int j = 0; j < ny; ++j) {  //overwrite prev_Vol_values, prev_Strike_grid
            K = yMin + j * dy;
            prev_Strike_grid[j] = K;
            prev_Vol_values[j] = mm_part(EPSILON, K); // INIZIALIZATION!! проверить ЧТО ОЧИСТИЛА PASS term ??? it is incide smiler EPSILON = expiry = 0
    }
    NumberFromFixLegSchedule = 2; // GET FROM PRODUCT for each swap 
    prev_s = 0.2; // FROM SMILER(expiry = ) get Swap value at  expiry = 0 !!!, term = xMax - T // todays swap spot swap   SPOT SWAP - TERM IS Fixed
    prev_dsdx = 0.01; // FROM SMILER(expiry = ) get dS/dx at expiry = 0, term = xMax - T// todays swap spot swap Use Market Curves 1/annnuity (....) 
                                //implement inside PRODUCT! or calculate here using market curves
    // -------------- END OF INITIALIZATION

    //--------MC containers
    std::vector<double> x_(NUMBER_OF_PARTICLES, 0);
    std::vector<double> y_(NUMBER_OF_PARTICLES, 0);
    std::vector<double> prev_s_(NUMBER_OF_PARTICLES, prev_s);
    std::vector<double> prev_dsdx_(NUMBER_OF_PARTICLES, prev_dsdx);
    std::vector<double> A_(NUMBER_OF_PARTICLES, 0);

    double random_val= 1;
    double sigma;
    double r;

    //---reserve for swap calculation
    double Tm;
    double tau_m;
    double P_tTm;
    double Level;
    double Weighted_Level;
    double f_tt;
    double f_tt_term;
    double P_tt_term;
    double sum;
    double THETA;

    for (int i = 1; i <= TIME_STEPS_PER_DELTA * nx; i++) {
        ti = xMin + i * particle_dt;
        term = xMax - ti;
        for(int j = 0; j < NUMBER_OF_PARTICLES; ++j) {

            double& x = x_[j]; //rename  particles correspond thier path ? or just write x_[i] isted of x
            double& y = y_[j];
            double& A = A_[j]; 
            //prev_Strike_grid, prev_Vol_values will be updated
            sigma = get_vol_from_smile(prev_Strike_grid, prev_Vol_values, prev_s_[j]) / prev_dsdx_[j];
            x = x + (y - LA * x) * particle_dt + sigma * std::sqrt(particle_dt) * random_val; //sqrt(dt)*W  state varible
            y = y + (sigma * sigma - 2 * LA * y)* particle_dt;  //state variable
            r = x + f0(ti); //noise, put in A = 
            A = A + r * particle_dt;
            //use NumberFromFixLegSchedule - is the number of payments over the term
            Level = 0;
            Weighted_Level = 0;
            tau_m = term/NumberFromFixLegSchedule;
            Tm = ti + tau_m;
            for(int m = 1; m < NumberFromFixLegSchedule + 1; m++){
                P_tTm = P0(Tm) / P0(ti) * std::exp(-G(ti, Tm) * x - 1/2 * G(ti, Tm) * G(ti, Tm) * y);
                Level += tau_m * P_tTm;
                Weighted_Level += tau_m * P_tTm * G(ti, Tm);
                Tm += tau_m;
            }

            //f_tt = f0(ti) + (x + G(ti, ti) * y); G(ti, ti) == 0 !!!
            f_tt_term = f0(ti+term) + std::exp(-LA*(term)) * (x + G(ti, ti+term) * y);
            P_tt_term = P0(ti+term) / P0(ti) * std::exp(-G(ti, ti+term) * x - 1/2 * G(ti, ti+term) * G(ti, ti+term) * y);

            prev_s_[j] = (1 - P_tt_term) / Level; //update
            prev_dsdx_[j] = -1.0/Level * (1 * G(ti, ti) - P_tt_term * G(ti, ti + term)) + prev_s_[j]/Level * Weighted_Level ; //update
        }

        if(i % TIME_STEPS_PER_DELTA == 0) {
            T = ti;
            yMin = 0.1; //Get from CDF using using smiler 
            yMax = 1; //Get from CDF
            dy = (yMax - yMin) / ny;  //not including right limit

            for (int l = 0; l < ny; ++l) {  //overwrite prev_Vol_values, prev_Strike_grid. // swapp lopps here?
                K = yMin + l * dy;
                
                sum = 0;

                for(int k = 0; k< NUMBER_OF_PARTICLES; k++) { //read x_ y_ ONLY!!

                    r = x_[k] + f0(ti);
                    Level = 0;
                    Weighted_Level = 0;
                    tau_m = term/NumberFromFixLegSchedule;
                    Tm = ti + tau_m; // !!
                    for(int m = 1; m < NumberFromFixLegSchedule + 1; m++){
                        P_tTm = P0(Tm) / P0(ti) * std::exp(-G(ti, Tm) * x_[k] -  0.5 * G(ti, Tm) * G(ti, Tm) * y_[k]);
                        Level += tau_m * P_tTm;
                        Weighted_Level += tau_m * P_tTm * G(ti, Tm);
                        Tm += tau_m;
                    }
                    f_tt = f0(ti) + (x_[k] + G(ti, ti) * y_[k]);  // G(ti, ti)  == 0 !!!
                    f_tt_term = f0(ti+term) + std::exp(-LA*(term)) * (x_[k] + G(ti, ti+term) * y_[k]);
                    P_tt_term = P0(ti+term) / P0(ti) * std::exp(-G(ti, ti+term) * x_[k] - 1/2 * G(ti, ti+term) * G(ti, ti+term) * y_[k]);
    
                    sum = sum + std::exp(-A_[k]) * (prev_s_[k] > K ? (f_tt - f_tt_term*P_tt_term - prev_s_[k]*(1 - P_tt_term)) : 0);
    

                }
                THETA = sum / NUMBER_OF_PARTICLES ;
                std::size_t raw_index = (int)(i / TIME_STEPS_PER_DELTA - 1);

                //----calculate mm_part here! and /CKK
                double dC_dKK = 1; /// WILL modify
                double volaa_in_2 = mm_part(T, K) + 2.0*THETA/dC_dKK; 
                double sigma = std::sqrt(std::max(0.0, volaa_in_2));
                vol_surface[raw_index][l] = sigma; 

                prev_Strike_grid[l] = K;
                prev_Vol_values[l] = vol_surface[raw_index][l] ; // INIZIALIZATION!! проверить ЧТО ОЧИСТИЛА PASS term ??? it is incide smiler EPSILON = expiry = 0
                grid[raw_index][l] = {T, K};
            
            }

          
        }

    }    //print grid and vol surface
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            std::cout << "T: " << grid[i][j].first << "K:" << grid[i][j].second << ", Vol: " << vol_surface[i][j] << std::endl;
        }
        std::cout << "-----------------------------" << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution time: " << elapsed.count() << " seconds\n";

    return 0;
}
