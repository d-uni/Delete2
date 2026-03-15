#include <iostream>
#include <string>
#include <utility>  // For std::move
#include <vector>
#include <algorithm>
#include <chrono>

////TRY WHEN WE FIX TERM === WRONG< WE want underlying cotermin. so we fix the final xMax, and simulate s(t, xMax-t) a,d not s(t, term)

int NUMBER_OF_PARTICLES = 100; //param by default = ...
int TIME_STEPS = 365; //param by default = ...
float EPSILON = 1e-5; //0.00001 put inside calibrater-private: static constexpr double EPSILON = 1e-9;

float LA = 1;


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

// we will change prev_Vol_values inside , so it is  NOT const 
std::vector<double> mc_part(double T, const std::vector<double>& Strike_grid, const std::vector<double>& prev_Vol_values, double prev_s, double prev_dsdx,  double term, int total_time_steps, int NumberFromFixLegSchedule) { //calculate THETA 
    double dt = T / total_time_steps;
    double random_val= 1;

    auto G = [](double t1, double t2)
    {
        return (1 - std::exp(-LA*(t2 - t1))) / LA;
    };

 /* double x = 0;
    double y = 0;
    double A = 0; */

    double ti;
    double THETA;
    double sum;
    double f_tt;
    double f_tt_term;
    double P_tt_term;
    double sigma;

    double r;
//---reserve for swap calculation
    double Tm;
    double tau_m;
    double P_tTm;
    double Level;
    double Weighted_Level;

    std::vector<double> x_(NUMBER_OF_PARTICLES, 0);
    std::vector<double> y_(NUMBER_OF_PARTICLES, 0);
    std::vector<double> A_(NUMBER_OF_PARTICLES, 0);
    std::vector<double> new_Vol_values(prev_Vol_values.size());
    
    for (int k = 0; k < Strike_grid.size(); k++) {
        for (int i = 1; i <= total_time_steps; ++i) { // NOT TAKING INTO ACCOUNT THAT t0 = EPSILON =! 0
            ti = i*dt;

            sum  = 0; 
            //все максимально убрать из цикла чтоь никаакие функции не вызывались 
            for(int j = 0; j < NUMBER_OF_PARTICLES; ++j) {
               double& x = x_[j]; //rename si particles correspond thier path ? or just write x_[i] isted of x
               double& y = y_[j];
               double& A = A_[j];
                sigma = get_vol_from_smile(Strike_grid, prev_Vol_values, prev_s) / prev_dsdx;
                x = x + (y - LA * x) * dt + sigma * random_val; //sqrt(dt)*W
                y = y + (sigma * sigma - 2 * LA * y)* dt;
                r = x + f0(ti);
                A = A + r*dt;

            //use NumberFromFixLegSchedule - is the number of payments over the term
                Level = 0;
                Weighted_Level = 0;
                tau_m = term/NumberFromFixLegSchedule;
                Tm = ti;
                for(int m = 1; m < NumberFromFixLegSchedule + 1; m++){
                    P_tTm = P0(Tm) / P0(ti) * std::exp(-G(ti, Tm) * x - 1/2 * G(ti, Tm) * G(ti, Tm) * y);
                    Level += tau_m * P_tTm;
                    Weighted_Level += tau_m * P_tTm * G(ti, Tm);
                    Tm += tau_m;
                }

                f_tt = f0(ti) + (x + G(ti, ti) * y);
                f_tt_term = f0(ti+term) + std::exp(-LA*(term)) * (x + G(ti, ti+term) * y);
                P_tt_term = P0(ti+term) / P0(ti) * std::exp(-G(ti, ti+term) * x - 1/2 * G(ti, ti+term) * G(ti, ti+term) * y);

                prev_s = (1 - P_tt_term) / Level; //update
                prev_dsdx = -1/Level * (1 * G(ti, ti) - P_tt_term) + prev_s/Level * Weighted_Level ; //update

                f_tt = f0(ti) + (x + G(ti, ti) * y);
                f_tt_term = f0(ti+term) + std::exp(-LA*(term)) * (x + G(ti, ti+term) * y);
                P_tt_term = P0(ti+term) / P0(ti) * std::exp(-G(ti, ti+term) * x - 1/2 * G(ti, ti+term) * G(ti, ti+term) * y);

                sum += std::exp(-A) * (prev_s > Strike_grid[k] ? (f_tt - f_tt_term*P_tt_term - prev_s*(1 - P_tt_term)) : 0);
                
            }
        
        }
        THETA = sum / NUMBER_OF_PARTICLES ;
        new_Vol_values[k] = THETA + mm_part(ti, Strike_grid[k]);

    }
    return new_Vol_values; //return vector of volatilities for each strike for given T, which is already updated with THETA
}


int main() {
    auto start = std::chrono::high_resolution_clock::now();
    //initialization (allocate memory once))
    double xMax = 10;
    double nx = 10;
    double xMin = 0;

    double yMax;  //Get from CDF
    double ny = 50;
    double yMin;  //Get from CDF 

    double dx = (xMax - xMin) / nx; //OR (xMax - xMin) / (nx - 1)
    double dy;

    //use statically allocated vectors
    std::vector<std::vector<std::pair<double,double>>> grid(nx, std::vector<std::pair<double,double>>(ny));
    std::vector<std::vector<double>> vol_surface(nx, std::vector<double>(ny));
    std::vector<double> Strike_grid(ny);
    std::vector<double> prev_Vol_values(ny);

    double T;  // build sigma(T, K) surface
    double K;
    double term;
    double particle_dt;
    double prev_s;
    double prev_dsdx ;
    int NumberFromFixLegSchedule;
    double mm_part_;
   // std::vector<double> mc_part_(ny);

    for (int i = 1; i <= nx; ++i) { 
        T = xMin + i * dx;                    // T = EXPIRY (of option ?) (xMax - T) = term 
        term = xMax - T;               //term of the option (time to maturity)

        //DETERMIN  yMax, yMin
        yMin = 0.1; //Get from CDF
        yMax = 1; //Get from CDF
        dy = (yMax - yMin) / ny;  //not including right limit

        // -------------- INITIALIZE  for MC part (same for any strike since we initialize the whole smile) + grid construction of STRIKE 
        //since we are alredy looping over T.

        for (int j = 0; j < ny; ++j) {  //overwrite prev_Vol_values
            K = yMin + j * dy;
            Strike_grid[j] = K;
            prev_Vol_values[j] = mm_part(EPSILON, K); // INIZIALIZATION!! проверить ЧТО ОЧИСТИЛА PASS term ??? it is incide smiler EPSILON = expiry = 0
        }
        NumberFromFixLegSchedule = 2; // GET FROM PRODUCT for each swap 
        prev_s = 0.2; // FROM SMILER(expiry = ) get Swap value at  expiry = 0 !!!, term = xMax - T // todays swap spot swap   SPOT SWAP - TERM IS Fixed
        prev_dsdx = 0.01; // FROM SMILER(expiry = ) get dS/dx at expiry = 0, term = xMax - T// todays swap spot swap Use Market Curves 1/annnuity (....) 
                                //implement inside PRODUCT! or calculate here using market curves
        // -------------- END OF INITIALIZATION


        //--------------Get THETA curve for all strikes 
        const std::vector<double>& mc_part_ = mc_part(T, Strike_grid, prev_Vol_values, prev_s, prev_dsdx, term, i * TIME_STEPS, NumberFromFixLegSchedule); //Preserve term only, ,model rolloing (spot) swap?
        
        for (int j = 0; j < ny; ++j) { 
            K = yMin + j * dy;
            mm_part_ = mm_part(T, K);
            vol_surface[i-1][j] = mm_part_ + mc_part_[j]; 
            grid[i-1][j] = {T, K};
        }
    }


    //we can put second loop  vol_surface[i-1][j] = mm_part_ + mc_part_[j];  insside mc_part_ and compute mm_part(T, K); in Strike Grid construction

/*
    //print grid and vol surface
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            std::cout << "T: " << grid[i][j].first << ", K:" << grid[i][j].second << ", Vol: " << vol_surface[i][j] << std::endl;
        }
        std::cout << "-----------------------------" << std::endl;
    }*/
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution time: " << elapsed.count() << " seconds\n";

    return 0;
}
