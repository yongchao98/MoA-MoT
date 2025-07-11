import pandas as pd
import numpy as np
import statsmodels.api as sm

def run_simulation():
    """
    Simulates data for three cases of treatment effect estimation
    and prints the results of OLS regressions.
    """
    # Set up a large population
    N = 100000
    np.random.seed(42)

    # X is pre-program income. Let's model it with a Lognormal distribution.
    X = np.random.lognormal(mean=9.5, sigma=0.8, size=N)
    
    # Potential outcome without treatment (Y0). It depends on pre-program income X.
    # The error term represents other unobserved factors.
    Y0 = 1000 + 0.9 * X + np.random.normal(loc=0, scale=1000, size=N)
    
    # The treatment effect is positive for everyone. Let's set it to a constant of 2000.
    treatment_effect = 2000
    
    # Potential outcome with treatment (Y1)
    Y1 = Y0 + treatment_effect

    print("--- Simulation Setup ---")
    print(f"Population size: {N}")
    print(f"True Average Treatment Effect: {treatment_effect}")
    print("Crucial Assumption: The treatment effect is positive for every individual.\n")

    # --- Case 1: D is randomly assigned ---
    print("--- Case 1: D is randomly assigned (simple randomization) ---")
    print("Model: Y ~ 1 + D")
    D1 = np.random.binomial(n=1, p=0.5, size=N)
    Y1_obs = Y0 * (1 - D1) + Y1 * (D1)
    
    df1 = pd.DataFrame({'Y': Y1_obs, 'D': D1, 'X': X})
    
    # Regression Y on a constant and D
    X_sm1 = sm.add_constant(df1[['D']])
    results1 = sm.OLS(df1['Y'], X_sm1).fit()
    
    b0 = results1.params['const']
    b1 = results1.params['D']
    
    print("Result: The coefficient on D represents the Average Treatment Effect (ATE).")
    print(f"Final Equation: Y = {b0:.2f} + {b1:.2f} * D")
    print(f"The coefficient on D is {b1:.2f}, which is positive, as expected.\n")
    

    # --- Case 2 & 3: D is assigned conditional on X ---
    # People with lower pre-program income (X) are more likely to be in the program.
    # This creates a correlation between D and X.
    propensity = np.where(X < np.median(X), 0.75, 0.25)
    D2 = np.random.binomial(n=1, p=propensity, size=N)
    Y2_obs = Y0 * (1 - D2) + Y1 * (D2)
    
    df2 = pd.DataFrame({'Y': Y2_obs, 'D': D2, 'X': X})
    
    # --- Case 2: Regress Y on D (omitting X) ---
    print("--- Case 2: D is assigned conditional on X (omitting X from regression) ---")
    print("Model: Y ~ 1 + D")
    
    # Regression Y on a constant and D
    X_sm2 = sm.add_constant(df2[['D']])
    results2 = sm.OLS(df2['Y'], X_sm2).fit()

    g0 = results2.params['const']
    g1 = results2.params['D']

    print("Result: This model suffers from omitted variable bias. The negative coefficient arises because the program targets lower-income individuals, and this pre-existing disadvantage is wrongly attributed to the program itself.")
    print(f"Final Equation: Y = {g0:.2f} + {g1:.2f} * D")
    print(f"The coefficient on D is {g1:.2f}. Even though the true effect is +{treatment_effect}, the coefficient is negative.\n")


    # --- Case 3: Regress Y on D and X ---
    print("--- Case 3: D is assigned conditional on X (controlling for X) ---")
    print("Model: Y ~ 1 + D + X")
    
    # Regression Y on a constant, D, and X
    X_sm3 = sm.add_constant(df2[['D', 'X']])
    results3 = sm.OLS(df2['Y'], X_sm3).fit()

    d0 = results3.params['const']
    d1 = results3.params['D']
    d2 = results3.params['X']

    print("Result: By controlling for X, we remove the omitted variable bias and recover an estimate of the treatment effect.")
    print(f"Final Equation: Y = {d0:.2f} + {d1:.2f} * D + {d2:.2f} * X")
    print(f"The coefficient on D is {d1:.2f}, which is positive and close to the true effect, as expected.\n")

if __name__ == '__main__':
    run_simulation()