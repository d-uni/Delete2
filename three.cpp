#include <iostream>
#include <string>
#include <utility>  // For std::move
#include <vector>
#include <algorithm>


int NUMBER_OF_PARTICLES = 1000; //param by default = ...
int TIME_STEPS = 1000; //param by default = ...
float EPSILON = 1e-5; //0.00001
float LA = 1;


double f0(double T) {
    return 0.12; // Market instantaneous FR
}

double P0(double T) {
    return std::exp(T * 0.12); // Market ZCB
}

double mm_part(double T, double K) {
    return 1; // here i have the correct market price
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

// we will change prev_Strike_grid and prev_Vol_values inside , so they are NOT const 
double mc_part(double T, double K, std::vector<double>& prev_Strike_grid, std::vector<double>&prev_Vol_values, double prev_s, double prev_dsdx,  double term, int total_time_steps, int NumberFromFixLegSchedule) { //calculate THETA 
    double dt = T / total_time_steps;
    double random_val= 1;

    auto G = [](double t1, double t2)
    {
        return (1 - std::exp(-LA*(t2 - t1))) / LA;
    };

    double x = 0;
    double y = 0;
    double A = 0;

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
    
    for (int i = 1; i <= total_time_steps; ++i) { // NOT TAKING INTO ACCOUNT THAT t0 = EPSILON =! 0
        ti = i*dt;
        for (int k = 0; k < prev_Strike_grid.size(); k++) {
            sum = 0;
            for(int j = 0; j < NUMBER_OF_PARTICLES; ++j) {
                sigma = get_vol_from_smile(prev_Strike_grid, prev_Vol_values, prev_s) / prev_dsdx;
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

                sum = std::exp(-A) * (prev_s > prev_Strike_grid[k] ? (f_tt - f_tt_term*P_tt_term - prev_s*(1 - P_tt_term)) : 0);
            
            }
            THETA = sum / NUMBER_OF_PARTICLES ;
            prev_Vol_values[k] = THETA + mm_part(ti, prev_Strike_grid[k]);
        }
    }

}


int main() {
    int xMax = 10; // NO INT
    int nx = 110;
    int xMin = 0;

    int yMax = 7;
    int ny = 10;
    int yMin = 1;

    double dx = (xMax - xMin) / (nx - 1);
    double dy;

    std::vector<std::vector<std::pair<double, double>>> grid; 
    // std::vector<std::vector<std::pair<double,double>>> grid(nx, std::vector<std::pair<double,double>>(ny));
    std::vector<std::vector<double>> vol_surface;  //    std::vector<std::vector<double>> vol_surface(nx, std::vector<double>(ny));
    std::vector<std::pair<double, double>> temp_vec;
    std::vector<double> temp_vec2;

    std::vector<double> prev_Strike_grid; //std::vector<double> prev_Strike_grid(ny);
    std::vector<double> prev_Vol_values; //std::vector<double> prev_Vol_values(ny);
    double T;
    double K;

    for (int i = 1; i <= nx; ++i) { //!!!!!!!!!!!!write <= if use static ellocation
        T = xMin + i * dx; // T - EXPIRY (of option ?) (xMax - T) - term 
        double term = xMax - T; //term of the option (time to maturity)

        temp_vec.clear(); // !! проверить ЧТО ОЧИСТИЛА 
        temp_vec2.clear(); // !! проверить ЧТО ОЧИСТИЛА


        // -------------- INITIALIZE  for MC part (same for any strike since we initialize the whole smile)
        prev_Strike_grid.clear(); // !! проверить ЧТО ОЧИСТИЛА
        prev_Vol_values.clear(); // !! проверить ЧТО ОЧИСТИЛА
    
        double particle_dt = T / (i * TIME_STEPS); //put it stratinto function write  // CHANGE HERE IF T/(TIME_STEPS)
        //write grid construction of STRIKE для Т = particle_dt один на всех, the maturity = delta  what is yMIN and yMAX?
        dy = (yMax - yMin) / (ny - 1);
        for (int j = 0; j < ny; ++j) {
            K = yMin + j * dy;
            prev_Strike_grid.push_back(K);
            prev_Vol_values.push_back(mm_part(EPSILON, K)); // !! проверить ЧТО ОЧИСТИЛА 
        }
        int NumberFromFixLegSchedule = 2; // GET FROM PRODUCT for each swap 
        double prev_s = 0.2; // FROM SMILER(expiry = ) get Swap value at  expiry = 0 !!!, term = xMax - T // todays swap spot swap   SPOT SWAP
        double prev_dsdx = 0.01; // FROM SMILER(expiry = ) get dS/dx at expiry = 0, term = xMax - T// todays swap spot swap Use Market Curves 1/annnuity (....) implement inside PRODUCT!
        // -------------- END OF INITIALIZATION
        

        //DETERMIN  yMax, yMin
        //double term = xMax - T; //term of the option (time to maturity)  USE FORWARD SWAP
        dy = (yMax - yMin) / (ny - 1);
        for (int j = 0; j < ny; ++j) {
            K = yMin + j * dy;
            temp_vec.push_back({T, K}); //grid[i][j] = {T, K}
                                        ////vol_surface[i][j] = mm_part_ + mc_part_
            
            double mm_part_ = mm_part(T, K);
            double mc_part_ = mc_part(T, K, prev_Strike_grid, prev_Vol_values, prev_s, prev_dsdx, term, i * TIME_STEPS, NumberFromFixLegSchedule); //Preserve term only, ,model rolloing (spot) swap!
            temp_vec2.push_back(mm_part_ + mc_part_);



            // How to make another name for function mc_part with fixed K and make it function of T only?
            //auto mc_part_T = [K](double T) { return mc_part(T, K); };
            //double mc_part_ = mc_part_T(T);
        }
        grid.push_back(temp_vec); 
        vol_surface.push_back(temp_vec2); 
    }

    return 0;
}
