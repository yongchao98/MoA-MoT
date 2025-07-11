import numpy as np
import pandas as pd
import statsmodels.api as sm

def run_simulation():
    """
    Simulates data to test the sign of the treatment effect coefficient
    in three different regression scenarios.
    """
    # Set parameters for the simulation
    np.random.seed(42)  # for reproducibility
    n_obs = 100000       # number of individuals
    true_treatment_effect = 2000  # a positive effect for everyone

    # 1. Generate base data (pre-program income X and potential outcomes)
    # X = pre-program income, from $20k to $70k
    X = np.random.uniform(20000, 70000, n_obs)
    
    # Potential outcome without treatment (Y0) depends on X + some random noise
    Y0 = 5000 + 1.1 * X + np.random.normal(0, 5000, n_obs)
    
    # Potential outcome with treatment (Y1) is Y0 + the true treatment effect
    # This ensures the treatment effect is positive for every single person.
    Y1 = Y0 + true_treatment_effect
    
    # --- Case 1: D is randomly assigned ---
    print("--- Case 1: Y on D (D is purely random) ---")
    # D is assigned with a 50% probability, independent of X
    D1 = np.random.binomial(1, 0.5, n_obs)
    # Observed income Y is Y1 if treated (D=1), Y0 if not (D=0)
    Y1_obs = Y1 * D1 + Y0 * (1 - D1)
    
    # Regression: Y ~ 1 + D
    exog1 = sm.add_constant(D1)
    model1 = sm.OLS(Y1_obs, exog1).fit()
    c1, d1_coef = model1.params
    
    print(f"The estimated coefficient for D is {d1_coef:.2f}, which is positive as expected.")
    print(f"Estimated Equation: Y = {c1:.2f} + {d1_coef:.2f} * D\n")

    # --- Case 2: D depends on X, but X is omitted from regression ---
    print("--- Case 2: Y on D (D depends on X, X omitted) ---")
    # D is assigned based on X: lower X means higher probability of treatment
    # This makes X a confounder.
    prob_d_depends_on_x = 1 / (1 + np.exp(0.0002 * (X - np.mean(X))))
    D2 = np.random.binomial(1, prob_d_depends_on_x, n_obs)
    
    # Observed income Y
    Y2_obs = Y1 * D2 + Y0 * (1 - D2)

    # Regression: Y ~ 1 + D
    exog2 = sm.add_constant(D2)
    model2 = sm.OLS(Y2_obs, exog2).fit()
    c2, d2_coef = model2.params
    
    print(f"The estimated coefficient for D is {d2_coef:.2f}.")
    if d2_coef < 0:
        print("It is NEGATIVE, demonstrating how omitted variable bias can reverse the sign.")
    print(f"Estimated Equation: Y = {c2:.2f} + {d2_coef:.2f} * D\n")
    
    # --- Case 3: D depends on X, and X is included in regression ---
    print("--- Case 3: Y on D and X (D depends on X, X included) ---")
    # Data is the same as in Case 2
    
    # Regression: Y ~ 1 + D + X
    df = pd.DataFrame({'Y': Y2_obs, 'D': D2, 'X': X})
    exog3 = sm.add_constant(df[['D', 'X']])
    model3 = sm.OLS(df['Y'], exog3).fit()
    c3, d3_coef, x3_coef = model3.params
    
    print(f"The estimated coefficient for D is {d3_coef:.2f}.")
    if d3_coef > 0:
        print("It is positive, showing that controlling for the confounder X recovers the true effect's sign.")
    print(f"Estimated Equation: Y = {c3:.2f} + {d3_coef:.2f} * D + {x3_coef:.2f} * X\n")

    print("Conclusion: The coefficient on D must be positive in Case 1 and Case 3.")

run_simulation()