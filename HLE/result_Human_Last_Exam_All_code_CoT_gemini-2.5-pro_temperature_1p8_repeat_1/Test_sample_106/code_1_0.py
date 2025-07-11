import numpy as np
import statsmodels.api as sm

def run_simulation():
    """
    Simulates data to demonstrate the effect of different regression models
    on estimating a treatment effect.
    """
    # --- 0. Set up parameters ---
    np.random.seed(42)  # for reproducibility
    n_individuals = 100000
    true_treatment_effect = 5  # Y(1) - Y(0) is always 5 (positive)
    
    print("This simulation demonstrates when the regression coefficient for a treatment 'D' must be positive.")
    print(f"We assume the true treatment effect on income 'Y' is always positive (+${true_treatment_effect}k for this simulation).\n")

    # --- 1. Generate underlying data ---
    # X = pre-treatment income. It's a strong predictor of post-treatment income.
    X = np.random.normal(loc=50, scale=10, size=n_individuals)
    # Potential outcomes: Y0 (income without treatment), Y1 (income with treatment)
    # They both depend on pre-treatment income X.
    noise = np.random.normal(loc=0, scale=5, size=n_individuals)
    Y0 = 20 + 0.8 * X + noise
    Y1 = Y0 + true_treatment_effect

    # --- 2. Case 1: D is randomly assigned ---
    print("--- Case 1: D is randomly assigned. We regress Y on D. ---")
    D_case1 = (np.random.rand(n_individuals) < 0.5).astype(int)
    Y_obs_case1 = Y1 * D_case1 + Y0 * (1 - D_case1)
    
    # Regression: Y ~ 1 + D
    reg_model_1 = sm.add_constant(D_case1)
    results_1 = sm.OLS(Y_obs_case1, reg_model_1).fit()
    const_1, coef_D_1 = results_1.params

    print(f"The resulting equation is: Y = {const_1:.2f} + {coef_D_1:.2f}*D")
    print(f"The coefficient on D is positive, as expected. It is very close to the true effect ({true_treatment_effect}).\n")

    # --- 3. Cases 2 & 3: D is assigned based on X (selection bias) ---
    print("--- Cases 2 & 3: D is assigned based on pre-treatment income X. ---")
    # People with lower income X are more likely to get treatment D=1.
    prob_d_equals_1 = 1 / (1 + np.exp((X - 50) / 5)) 
    D_cases_2_3 = (np.random.rand(n_individuals) < prob_d_equals_1).astype(int)
    Y_obs_cases_2_3 = Y1 * D_cases_2_3 + Y0 * (1 - D_cases_2_3)

    # --- 4. Case 2: Conditional assignment, but regression omits X ---
    print("\n--- Case 2: Regress Y on D (omitting X). ---")
    # Regression: Y ~ 1 + D
    reg_model_2 = sm.add_constant(D_cases_2_3)
    results_2 = sm.OLS(Y_obs_cases_2_3, reg_model_2).fit()
    const_2, coef_D_2 = results_2.params

    print(f"The resulting equation is: Y = {const_2:.2f} + {coef_D_2:.2f}*D")
    print("The coefficient on D is negative. This is due to omitted variable bias: the program appears harmful because it's mainly taken by people who would have had lower incomes anyway.\n")

    # --- 5. Case 3: Conditional assignment, regression includes X ---
    print("\n--- Case 3: Regress Y on D and X. ---")
    # Regression: Y ~ 1 + D + X
    reg_model_3 = sm.add_constant(np.column_stack((D_cases_2_3, X)))
    results_3 = sm.OLS(Y_obs_cases_2_3, reg_model_3).fit()
    const_3, coef_D_3, coef_X_3 = results_3.params

    print(f"The resulting equation is: Y = {const_3:.2f} + {coef_D_3:.2f}*D + {coef_X_3:.2f}*X")
    print("The coefficient on D is positive. By controlling for X, we remove the selection bias and correctly estimate the positive treatment effect.\n")

if __name__ == '__main__':
    run_simulation()
<<<E>>>