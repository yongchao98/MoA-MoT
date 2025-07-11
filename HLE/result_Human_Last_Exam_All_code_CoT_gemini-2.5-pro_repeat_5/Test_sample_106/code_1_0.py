import numpy as np
import statsmodels.api as sm

def run_simulation():
    """
    Simulates three scenarios to test the sign of a treatment effect coefficient.
    
    Case 1: D is fully randomized. Regression: Y ~ D.
    Case 2: D is conditionally randomized on X. Regression: Y ~ D.
    Case 3: D is conditionally randomized on X. Regression: Y ~ D + X.
    
    The true treatment effect is set to be positive for all individuals.
    """
    # 0. Set up the population and parameters
    np.random.seed(42)
    N = 100000  # Population size
    
    # Generate pre-treatment income X
    X = np.random.lognormal(mean=np.log(20000), sigma=0.5, size=N)
    
    # Define potential outcomes Y(0) and Y(1)
    # Y(0) is income without the program. It depends strongly on pre-treatment income X.
    # We add some noise to make it more realistic.
    Y0 = 1.5 * X + np.random.normal(0, 1000, size=N)
    
    # The treatment effect is a positive constant for everyone.
    # This ensures the core assumption (Y(1) > Y(0)) holds for all.
    treatment_effect = 2000
    Y1 = Y0 + treatment_effect
    
    print("--- Simulation Setup ---")
    print(f"True Treatment Effect for everyone: {treatment_effect}\n")

    # --- Case 1: D is randomly assigned ---
    # D is assigned with a 50% probability, independent of X
    D_case1 = np.random.binomial(1, 0.5, size=N)
    Y_case1 = D_case1 * Y1 + (1 - D_case1) * Y0
    
    # Regress Y on a constant and D
    X_case1_reg = sm.add_constant(D_case1)
    model1 = sm.OLS(Y_case1, X_case1_reg).fit()
    coef_D_case1 = model1.params[1]
    
    print("--- Case 1: Y ~ D (D is fully random) ---")
    print(f"The estimated coefficient on D is: {coef_D_case1:.2f}")
    print("Result: The coefficient is positive, as expected.\n")

    # --- Case 2 & 3: D is randomly assigned conditional on X ---
    # Probability of treatment is higher for lower X (program targets low-income)
    # The formula creates a smooth probability transition based on X
    propensity_score = 1 / (1 + np.exp(0.0002 * (X - np.median(X))))
    D_case23 = np.random.binomial(1, propensity_score, size=N)
    Y_case23 = D_case23 * Y1 + (1 - D_case23) * Y0
    
    # --- Case 2: Omitting the control variable X ---
    # Regress Y on a constant and D only
    X_case2_reg = sm.add_constant(D_case23)
    model2 = sm.OLS(Y_case23, X_case2_reg).fit()
    coef_D_case2 = model2.params[1]

    print("--- Case 2: Y ~ D (D depends on X, but X is omitted) ---")
    print(f"The estimated coefficient on D is: {coef_D_case2:.2f}")
    print("Result: The coefficient is negative due to omitted variable bias.\n")

    # --- Case 3: Including the control variable X ---
    # Regress Y on a constant, D, and X
    X_case3_reg = sm.add_constant(np.column_stack((D_case23, X)))
    model3 = sm.OLS(Y_case23, X_case3_reg).fit()
    coef_D_case3 = model3.params[1]

    print("--- Case 3: Y ~ D + X (D depends on X, and X is controlled for) ---")
    print(f"The estimated coefficient on D is: {coef_D_case3:.2f}")
    print("Result: The coefficient is positive and close to the true effect.\n")

run_simulation()
<<<E>>>