#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cassert>
#include "matrix.h"
#include "mex.h"


using namespace std;

int main(void);

void my_accumulate(vector<double>& values, vector<double>& output)
{
    output.resize(values.size());
    fill(output.begin(), output.end(), 0);
	
    int n = static_cast<int> (values.size());
	output[0] = values[0];
    for(int i = 1; i < n; i++)
    {
        output[i] = output[i - 1] + values[i];
    }
}

double select_best_delta(int n1, int n2, double beta_P1, double beta_P1_1, double beta_P2, double beta_P2_1, double sum_k_R1, double sum_k_R2, double lambda, double points_to_validate[4])
{
    double delta = 0.;
    double best_J = FLT_MAX;
    for(int i = 0; i < 4; i++){
        double x = points_to_validate[i];
        double next_J = x * (sum_k_R1 / n1 - sum_k_R2 / n2) + 
                       lambda * (((beta_P1 + x / n1) != (beta_P1_1 - x / n2)) + 
								 ((beta_P2 - x / n2) != (beta_P2_1 + x / n1)));

        if (next_J < best_J){
            best_J = next_J;
            delta = x;
        }
    }

    return delta;
}


void update_betas(vector<double>& betas, int P1, int P2, double delta, int n1, int n2)
{
    int i = 0;
	int n = n1 + n2;
    vector<double>::iterator it=betas.begin();
    for(; i < P2; i++){
        *it += delta / n1;
        ++it;
    }
    for(; i < P1; i++){
        *it -= delta / n2;
        ++it;
    }
    for(; i < n; i++){
        *it += delta / n1;
        ++it;
    }
}

double globalJ(vector<double>& curvature, vector<double>& betas, double lambda)
{
    int i = 0;
	int n = static_cast<int> (curvature.size());
    double accum_curvature = 0.;
    double accum_regularization = 0.;

    for (int i = 0; i < n; i++){
        accum_curvature += curvature[i] * betas[i];
        accum_regularization += (betas[i] != betas[n - 1 - (n - i) * ((i - 1) >= 0)]);
    }

    return accum_curvature + lambda * accum_regularization;
}

mxArray * getMexArray(vector<double>& v)
{
	mxArray * mx = mxCreateDoubleMatrix(1, v.size(), mxREAL);
	copy(v.begin(), v.end(), mxGetPr(mx));
	return mx;
}


void optimizer_l0_m(double lambda,
					double max_value,
					int max_iterations,
					double * k_ptr,
					double * b_ptr,
					vector<double>& final_betas)
{

	int n = (int) final_betas.size();

	vector<double> curvature(k_ptr, k_ptr + n);
	vector<double> betas(b_ptr, b_ptr + n);

	vector<double> cumulative_k;
	double points_to_validate[4];

    double mean_value = static_cast<double> (accumulate(betas.begin(), betas.end(), 0.0));
    mean_value /= n;

    my_accumulate(curvature, cumulative_k);
    double total_k = static_cast<double> (accumulate(curvature.begin(), curvature.end(), 0.0));

	// initial solution loss
    double weighted_k = 0;
	double reg_term = 0;

	for(int i = 0; i < n; i++){
		weighted_k += curvature[i] * betas[i];
		reg_term += (betas[i] != betas[(i + 1) % n]);
	}
	double bestJ = weighted_k + lambda * reg_term;
	vector<double> base_betas = betas;

	// for each round of optimization
    for(int round = 0; round < max_iterations; round++){

		double this_it_bestJ = FLT_MAX;
		vector<double> best_betas(n);

        pair<pair<int, int>, double> best_update;

		// for each pair of points
        for(int P1 = 1; P1 < n; P1++){
            for (int P2 = 0; P2 < P1; P2++){

                int n1 = n - P1 + P2;
				int n2 = n - n1;				

				double sum_k_R2 = cumulative_k[P1] - cumulative_k[P2];
				double sum_k_R1 = total_k - sum_k_R2;


				// Finding the minimum and maximum gamma values admissible for the current pair of points

                double min_b1=FLT_MAX, max_b1=FLT_MIN, 
					  min_b2=FLT_MAX, max_b2=FLT_MIN;

                int i = 0;
                vector<double>::iterator it=base_betas.begin();
                for(; i < P2; i++){
                    min_b1 = min(min_b1, *it);
                    max_b1 = max(max_b1, *it);
                    ++it;
                }
                for(; i < P1; i++){
                    min_b2 = min(min_b2, *it);
                    max_b2 = max(max_b2, *it);
                    ++it;
                }
                for(; i < n; i++){
                    min_b1 = min(min_b1, *it);
                    max_b1 = max(max_b1, *it);
                    ++it;
                }

				
				double min_delta = max(-n1 * min_b1, n2 * (max_b2 - max_value));
				double max_delta = min(n1 * (max_value - max_b1), n2 * min_b2);

				// one of the following values is the optimum gamma
                points_to_validate[2] = min_delta;
                points_to_validate[3] = max_delta;
                points_to_validate[1] = ((n1*n2) / (n1 + n2)) * (base_betas[P1 - 1] - base_betas[P1]);
                points_to_validate[0] = ((n1*n2) / (n1 + n2)) * (base_betas[P2] - base_betas[n - 1 - (n - P2) * ((P2-1)>=0)]);

				double delta = select_best_delta(n1,
												n2, 
												base_betas[P1], 
												base_betas[P1-1],
												base_betas[P2],
												base_betas[n - 1 - (n - P2) * ((P2-1)>=0)],
												sum_k_R1,
												sum_k_R2,
												lambda,
												points_to_validate);

				if (delta < min_delta)
					delta = min_delta;

				if (delta > max_delta)
					delta = max_delta;

				vector<double> this_config_betas = base_betas;
                update_betas(this_config_betas, P1, P2, delta, n1, n2);
				double this_configJ = globalJ(curvature, this_config_betas, lambda);


                if (this_configJ < this_it_bestJ){
                    //best_update = make_pair(make_pair(P1, P2),
                    //                        delta);

					this_it_bestJ = this_configJ;
					best_betas = this_config_betas;
                }

			}
		}

        //int P1_s = best_update.first.first;
        //int P2_s = best_update.first.second;
        //int n2_s = (P1_s - P2_s);
        //int n1_s = n - n2_s;
        //double delta = best_update.second;

        if (this_it_bestJ < bestJ){
            bestJ = this_it_bestJ;
        }
		else
			break;

		base_betas = best_betas;

	}

	final_betas = base_betas;

    return;
}


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	if (nrhs != 5) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", "MEXCPP requires five input arguments.");
	} else if (nlhs > 1) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout", "MEXCPP requires a single output argument.");
	}

	/* Input declaration and assignment*/
	double lambda = (double) mxGetScalar(prhs[0]);
	double max_value = (double) mxGetScalar(prhs[1]);
	int max_iterations = (int) mxGetScalar(prhs[2]);
	double * k_ptr = mxGetPr(prhs[3]);
	double * b_ptr = mxGetPr(prhs[4]);


	/* Output declaration and assignment*/
	int array_size = (int) mxGetN(prhs[3]);
	vector<double> final_betas(array_size);

	optimizer_l0_m(lambda, max_value, max_iterations, k_ptr, b_ptr, final_betas);

	plhs[0] = getMexArray(final_betas);

	return;
}
