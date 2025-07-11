import numpy as np
from scipy.stats import chi2

def perform_hansen_j_test():
    """
    This function simulates data for an instrumental variable model,
    performs a 2SLS estimation, and then runs a Hansen J-test (Sargan variant)
    to check the validity of the instruments.

    Scenario:
    - Y: outcome
    - D: endogenous variable
    - Z1, Z2: Valid instruments
    - Z3: Invalid instrument (correlated with the structural error term)
    
    The J-test tests the null hypothesis that all instruments are valid.
    Since one instrument (Z3) is invalid, we expect the test to reject the null.
    """
    # 1. Simulate Data
    np.random.seed(42)
    n = 1000  # Sample size
    
    # Instruments
    Z1 = np.random.normal(0, 1, n)
    Z2 = np.random.normal(0, 1, n)
    Z3 = np.random.normal(0, 1, n)

    # Structural error term 'u'. We make it correlated with Z3 to make Z3 invalid.
    u = 0.8 * Z3 + np.random.normal(0, 1, n)
    
    # First stage equation for the endogenous variable D
    # D = f(Z1, Z2, Z3, v)
    v = np.random.normal(0, 1, n)
    # D is made continuous for simplicity in the 2SLS calculation, but the principle is the same.
    D = 1.5 * Z1 + 1.5 * Z2 + 0.5 * Z3 + v

    # Structural equation (second stage)
    # Y = beta_0 + beta_1 * D + u
    true_beta_1 = 2.0
    Y = 5 + true_beta_1 * D + u
    
    print("--- Data Simulation ---")
    print(f"Simulated data with {n} observations.")
    print("Instruments Z1 and Z2 are valid.")
    print("Instrument Z3 is invalid (correlated with the structural error term 'u').\n")

    # 2. Perform 2SLS and the J-Test
    # We will test the validity of {Z1, Z2, Z3} jointly.
    # If the test rejects, and we maintain that Z1 is valid,
    # the invalidity must come from Z2 or Z3.

    # Setup matrices
    # Endogenous regressors matrix (intercept and D)
    X = np.vstack([np.ones(n), D]).T 
    # Instrument matrix (intercept, Z1, Z2, Z3)
    Z = np.vstack([np.ones(n), Z1, Z2, Z3]).T 

    # --- Manual 2SLS ---
    # First stage: Regress D on Z to get D_hat
    # Since we have D (the single endogenous variable), we focus on that column
    d_vec = X[:, 1]
    d_hat_params = np.linalg.lstsq(Z, d_vec, rcond=None)[0]
    d_hat = Z @ d_hat_params
    
    # Create the matrix with predicted endogenous variable
    X_hat = np.vstack([np.ones(n), d_hat]).T

    # Second stage: Regress Y on X_hat to get 2SLS estimates
    beta_2sls = np.linalg.lstsq(X_hat, Y, rcond=None)[0]
    
    print("--- 2SLS Estimation ---")
    print(f"Estimated beta_0 (intercept): {beta_2sls[0]:.4f}")
    print(f"Estimated beta_1 (effect of D on Y): {beta_2sls[1]:.4f} (True value is {true_beta_1})")

    # 3. Calculate Hansen J-Statistic (Sargan's Test)
    # Get residuals from the original structural equation using 2SLS estimates
    residuals = Y - X @ beta_2sls

    # Regress residuals on the instruments Z
    res_on_z_params = np.linalg.lstsq(Z, residuals, rcond=None)[0]
    predicted_residuals = Z @ res_on_z_params
    
    # Calculate R-squared for this auxiliary regression
    ss_res = np.sum((residuals - predicted_residuals)**2)
    ss_tot = np.sum((residuals - np.mean(residuals))**2)
    r_squared = 1 - (ss_res / ss_tot)

    # J-statistic = n * R^2
    j_statistic = n * r_squared
    
    # Degrees of freedom = (number of instruments) - (number of endogenous variables)
    num_instruments = Z.shape[1] - 1 # Subtract one for the intercept in our matrix
    num_endog_vars = X.shape[1] - 1  # Subtract one for the intercept
    df = num_instruments - num_endog_vars
    
    # Calculate p-value from the chi-squared distribution
    p_value = 1 - chi2.cdf(j_statistic, df)
    
    print("\n--- Hansen J-Test for Overidentifying Restrictions ---")
    print("Null Hypothesis (H0): All instruments are valid (uncorrelated with errors).")
    print(f"Number of instruments: {num_instruments}")
    print(f"Number of endogenous variables: {num_endog_vars}")
    print(f"J-Statistic: {j_statistic:.4f}")
    print(f"Degrees of Freedom: {df}")
    print(f"P-value: {p_value:.4f}")

    # Conclusion
    alpha = 0.05
    if p_value < alpha:
        print(f"\nConclusion: Reject the null hypothesis at the {alpha*100}% significance level.")
        print("The test suggests that at least one instrument is invalid.")
    else:
        print(f"\nConclusion: Fail to reject the null hypothesis at the {alpha*100}% significance level.")
        print("The test does not provide evidence against the validity of the instruments.")

perform_hansen_j_test()
<<<B>>>