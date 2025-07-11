import numpy as np
import statsmodels.api as sm

def run_simulation_and_print_results():
    """
    Simulates data for the three cases and prints the regression results.
    """
    # Set seed for reproducibility
    np.random.seed(42)
    
    # 1. Setup a large population
    N = 200000
    
    # Generate X (pre-treatment income), ensuring it's positive
    # Use a log-normal distribution for a realistic income-like variable
    X = np.random.lognormal(mean=10.5, sigma=0.4, size=N)
    
    # Generate an unobserved error term
    u = np.random.normal(0, 5000, size=N)

    # 2. Define potential outcomes and a positive treatment effect (TE) for everyone
    # Let the TE depend on X, but ensure it's always positive
    # TE = Y(1) - Y(0)
    treatment_effect = 3000 + 0.05 * X 
    
    # Potential outcome if not treated (Y0)
    Y0 = 20000 + 1.2 * X + u
    
    # Potential outcome if treated (Y1)
    Y1 = Y0 + treatment_effect
    
    print("--- Data Generation Setup ---")
    print(f"Population size N = {N}")
    print(f"Average Treatment Effect (ATE) in population: {np.mean(treatment_effect):.2f}")
    print("The individual treatment effect is constructed to be positive for everyone.")
    print("-" * 30 + "\n")

    # --- Case 1: D is randomly assigned ---
    print("--- Case 1: D is randomly assigned. Regress Y on D. ---")
    # D is assigned with 50% probability
    D_case1 = (np.random.uniform(size=N) < 0.5).astype(int)
    # Observed outcome Y is a combination of potential outcomes based on treatment status
    Y_case1 = D_case1 * Y1 + (1 - D_case1) * Y0
    
    # Run regression Y ~ D
    X_case1_reg = sm.add_constant(D_case1)
    model_case1 = sm.OLS(Y_case1, X_case1_reg).fit()
    coef_d_case1 = model_case1.params[1]
    
    print("Regression: Y = b0 + b1*D")
    print(f"Estimated equation: Y = {model_case1.params[0]:.2f} + {coef_d_case1:.2f}*D")
    print(f"The coefficient on D is {coef_d_case1:.2f}, which is positive.")
    print("Conclusion: As expected, the coefficient approximates the ATE because D is randomized.\n")

    # --- Cases 2 & 3: D is randomly assigned conditional on X ---
    print("--- Cases 2 & 3: D assigned conditional on X. ---")
    # Let D be more likely for individuals with lower pre-treatment income (X)
    # This creates a selection bias.
    # We use a logistic function for the probability
    median_x = np.median(X)
    prob_d_equals_1 = 1 / (1 + np.exp(0.00005 * (X - median_x)))
    D_case23 = (np.random.uniform(size=N) < prob_d_equals_1).astype(int)
    
    # Observed outcome Y
    Y_case23 = D_case23 * Y1 + (1 - D_case23) * Y0

    # Verify that selection on X occurred
    mean_x_treated = np.mean(X[D_case23 == 1])
    mean_x_untreated = np.mean(X[D_case23 == 0])
    print(f"Mean X for treated group: {mean_x_treated:.2f}")
    print(f"Mean X for untreated group: {mean_x_untreated:.2f}")
    print("As designed, the treated group has a lower average pre-treatment income (X).\n")
    
    # --- Case 2: Regress Y on D (omitting X) ---
    print("--- Case 2: Regress Y on D (omitting X). ---")
    X_case2_reg = sm.add_constant(D_case23)
    model_case2 = sm.OLS(Y_case23, X_case2_reg).fit()
    coef_d_case2 = model_case2.params[1]
    
    print("Regression: Y = b0 + b1*D")
    print(f"Estimated equation: Y = {model_case2.params[0]:.2f} + {coef_d_case2:.2f}*D")
    print(f"The coefficient on D is {coef_d_case2:.2f}, which is negative.")
    print("Conclusion: The positive treatment effect is overwhelmed by the negative selection bias (treated people had much lower pre-treatment income X, which also leads to lower post-treatment income Y).\n")
    
    # --- Case 3: Regress Y on D and X ---
    print("--- Case 3: Regress Y on D and X. ---")
    X_case3_reg = sm.add_constant(np.column_stack((D_case23, X)))
    model_case3 = sm.OLS(Y_case23, X_case3_reg).fit()
    coef_d_case3 = model_case3.params[1]

    print("Regression: Y = b0 + b1*D + b2*X")
    print(f"Estimated equation: Y = {model_case3.params[0]:.2f} + {coef_d_case3:.2f}*D + {model_case3.params[2]:.2f}*X")
    print(f"The coefficient on D is {coef_d_case3:.2f}, which is positive.")
    print("Conclusion: After controlling for X, the coefficient on D correctly estimates the ATE and is positive.\n")

if __name__ == '__main__':
    run_simulation_and_print_results()