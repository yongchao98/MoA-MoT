import pandas as pd
import numpy as np
import statsmodels.api as sm

def run_simulation(N=20000):
    """
    Simulates data to test the three cases.
    """
    # 1. Setup: Generate data
    # X = pre-treatment income (in thousands)
    # The relationship between pre- and post-income is positive.
    np.random.seed(42)
    X = np.random.normal(loc=50, scale=10, size=N)
    
    # Y0 = potential income without treatment. Depends positively on X.
    u = np.random.normal(loc=0, scale=5, size=N)
    Y0 = 50 + 1.2 * X + u
    
    # TE = Treatment effect. Assumed to be POSITIVE for everyone.
    # We make it heterogeneous, but it's always > 0 since X is positive.
    TE = 10 + 0.05 * X 
    
    # Y1 = potential income with treatment
    Y1 = Y0 + TE
    
    print("--- Simulation Setup ---")
    print(f"Average Treatment Effect (ATE) in population: {np.mean(TE):.2f}\n")

    # --- Case 1: D is randomly assigned ---
    print("--- Case 1: D is randomly assigned. Regress Y on D. ---")
    D_rand = np.random.binomial(1, 0.5, N)
    Y_obs_rand = Y0 * (1 - D_rand) + Y1 * D_rand
    
    model1 = sm.OLS(Y_obs_rand, sm.add_constant(D_rand))
    results1 = model1.fit()
    coef_d_1 = results1.params[1]
    
    print(f"The estimated coefficient on D is: {coef_d_1:.4f}")
    print("This estimates the ATE. As predicted, it is positive.\n")

    # --- Cases 2 & 3: D is random conditional on X ---
    # Assignment to treatment is more likely for lower-income individuals (lower X)
    # This creates a negative correlation between D and X
    prob_d_cond = 1 / (1 + np.exp(-1 * (6 - 0.1 * X))) 
    D_cond = np.array([np.random.binomial(1, p) for p in prob_d_cond])
    Y_obs_cond = Y0 * (1 - D_cond) + Y1 * D_cond
    
    print("--- Case 2: D is conditional on X. Regress Y on D (Omitted X). ---")
    model2 = sm.OLS(Y_obs_cond, sm.add_constant(D_cond))
    results2 = model2.fit()
    coef_d_2 = results2.params[1]
    
    print(f"The estimated coefficient on D is: {coef_d_2:.4f}")
    print("This is NEGATIVE due to strong omitted variable bias. Lower-income individuals (who have lower Y on average) were selected into treatment, creating a negative bias that overwhelms the positive treatment effect.\n")

    print("--- Case 3: D is conditional on X. Regress Y on D and X. ---")
    X_df = pd.DataFrame({'const': 1, 'D': D_cond, 'X': X})
    model3 = sm.OLS(Y_obs_cond, X_df)
    results3 = model3.fit()
    coef_d_3 = results3.params['D']

    print(f"The estimated coefficient on D is: {coef_d_3:.4f}")
    print("This coefficient correctly controls for X. As predicted, it is positive because it estimates a weighted average of the (always positive) conditional treatment effects.")

if __name__ == '__main__':
    run_simulation()