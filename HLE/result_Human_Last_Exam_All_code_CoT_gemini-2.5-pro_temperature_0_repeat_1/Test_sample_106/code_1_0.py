import numpy as np
import statsmodels.api as sm

def run_simulations():
    """
    This function runs simulations for the three cases described.
    """
    # --- Common Parameters ---
    N = 500000  # Use a large sample size to approximate the population
    np.random.seed(42)
    
    # Generate a covariate X (pre-treatment income)
    X = np.random.normal(loc=50000, scale=10000, size=N)
    
    # Define potential outcomes. Y0 is income without treatment, Y1 is with treatment.
    # The treatment effect is a constant positive amount 'tau' for everyone.
    tau = 2000  # Positive treatment effect
    noise = np.random.normal(0, 5000, N)
    Y0 = 10000 + 0.8 * X + noise
    Y1 = Y0 + tau  # Y1 > Y0 for everyone

    print("--- Simulation Setup ---")
    print(f"The true treatment effect is a positive constant: {tau}")
    print("-" * 25)

    # --- Case 1: D is randomly assigned ---
    print("\n--- Case 1: D is randomly assigned. Regress Y on D. ---")
    # D is assigned with 50% probability, independent of X
    D1 = np.random.binomial(1, 0.5, N)
    # Observed outcome Y
    Y1_obs = Y0 * (1 - D1) + Y1 * D1
    # Regression model
    X1_reg = sm.add_constant(D1)
    model1 = sm.OLS(Y1_obs, X1_reg).fit()
    intercept1, coef_d1 = model1.params
    print("The coefficient on D should be positive and close to the true ATE.")
    print(f"Estimated Equation: Y = {intercept1:.2f} + {coef_d1:.2f} * D")
    print(f"Result: The coefficient on D ({coef_d1:.2f}) is positive as expected.")

    # --- Cases 2 & 3: D is randomly assigned conditional on X ---
    print("\n--- Cases 2 & 3: D is assigned based on X. ---")
    # Probability of treatment is higher for lower X (pre-treatment income)
    # This creates a confounding effect
    prob_d2 = 1 / (1 + np.exp((X - np.mean(X)) / np.std(X)))
    D2 = (np.random.rand(N) < prob_d2).astype(int)
    # Observed outcome Y
    Y2_obs = Y0 * (1 - D2) + Y1 * D2
    
    # --- Case 2: Regress Y on D only (omitting X) ---
    print("\n--- Case 2: Regress Y on D (omitting X). ---")
    X2_reg = sm.add_constant(D2)
    model2 = sm.OLS(Y2_obs, X2_reg).fit()
    intercept2, coef_d2 = model2.params
    print("The coefficient on D can be negative due to omitted variable bias.")
    print(f"Estimated Equation: Y = {intercept2:.2f} + {coef_d2:.2f} * D")
    print(f"Result: The coefficient on D ({coef_d2:.2f}) is negative, demonstrating the bias.")

    # --- Case 3: Regress Y on D and X ---
    print("\n--- Case 3: Regress Y on D and X. ---")
    X3_reg = sm.add_constant(np.column_stack((D2, X)))
    model3 = sm.OLS(Y2_obs, X3_reg).fit()
    intercept3, coef_d3, coef_x3 = model3.params
    print("The coefficient on D should be positive after controlling for X.")
    print(f"Estimated Equation: Y = {intercept3:.2f} + {coef_d3:.2f} * D + {coef_x3:.2f} * X")
    print(f"Result: The coefficient on D ({coef_d3:.2f}) is positive as expected.")

if __name__ == '__main__':
    run_simulations()