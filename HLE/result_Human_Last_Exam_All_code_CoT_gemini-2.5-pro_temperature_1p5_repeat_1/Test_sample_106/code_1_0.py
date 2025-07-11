import numpy as np
import pandas as pd
import statsmodels.api as sm

def run_simulation():
    """
    This simulation demonstrates the properties of the regression coefficient for a treatment
    effect under three different scenarios as described in the problem.
    """
    # Set up the simulation parameters for a large population
    np.random.seed(42) # for reproducibility
    n_population = 100000

    # --- 1. Define the Population and Potential Outcomes ---
    # X: pre-program income (in thousands). This is our control variable.
    X = np.random.normal(loc=50, scale=10, size=n_population)

    # Y(0): Potential outcome (post-program income) without treatment.
    # It depends on pre-program income X.
    error_y0 = np.random.normal(loc=0, scale=5, size=n_population)
    Y0 = 50 + 1.0 * X + error_y0

    # Treatment effect (tau): Positive for everyone, as per the problem statement.
    # We assume a constant treatment effect of 5 (i.e., $5,000).
    tau = 5.0

    # Y(1): Potential outcome with treatment.
    Y1 = Y0 + tau

    # Create a pandas DataFrame to hold our population data
    df = pd.DataFrame({'X': X, 'Y0': Y0, 'Y1': Y1})

    print("--- Simulation Setup ---")
    print(f"Population size: {n_population}")
    print(f"Treatment effect (tau) is constant and positive: {tau}")
    print("Potential outcome Y0 (no treatment) is positively correlated with X (pre-income).\n")

    # --- 2. Case 1: D is randomly assigned ---
    print("--- Case 1: D is randomly assigned ---")
    print("In this scenario, treatment D is independent of X.")
    D_case1 = np.random.binomial(1, 0.5, n_population)
    Y_case1 = D_case1 * df['Y1'] + (1 - D_case1) * df['Y0']

    # Regression: Y on a constant and D
    X_regr_case1 = sm.add_constant(D_case1)
    model_case1 = sm.OLS(Y_case1, X_regr_case1).fit()
    coef_D_case1 = model_case1.params[1]

    print("Running regression: Y ~ 1 + D")
    print(f"Estimated coefficient on D: {coef_D_case1:.4f}")
    print("Explanation: With random assignment, the coefficient on D estimates the Average Treatment Effect (ATE).")
    print("Since the treatment effect is positive for everyone, the ATE must be positive.")
    print("Conclusion for Case 1: The coefficient MUST be positive.\n")


    # --- 3. Setup for Case 2 & 3: D is assigned based on X ---
    print("--- Setup for Case 2 & 3: D is assigned based on X ---")
    # Here, treatment is targeted at individuals with lower pre-program income (X < 50).
    # This induces a correlation between D and X, leading to selection bias if X is omitted.
    print("Treatment is targeted at individuals with pre-income X < 50.")
    D_case23 = (df['X'] < 50).astype(int)
    Y_case23 = D_case23 * df['Y1'] + (1 - D_case23) * df['Y0']
    # Add these to a new DataFrame for clarity
    df_case23 = pd.DataFrame({'X': X, 'D': D_case23, 'Y': Y_case23, 'Y0': Y0})

    # --- 4. Case 2: D is assigned based on X, regress Y on D only ---
    print("\n--- Case 2: D is assigned based on X, regression Y ~ 1 + D ---")
    print("Here we omit the control variable X, creating omitted variable bias.")
    X_regr_case2 = sm.add_constant(df_case23['D'])
    model_case2 = sm.OLS(df_case23['Y'], X_regr_case2).fit()
    coef_D_case2 = model_case2.params[1]

    print("Running regression: Y ~ 1 + D")
    print(f"Estimated coefficient on D: {coef_D_case2:.4f}")
    # We can explicitly calculate the selection bias in our simulation
    treated_group_y0_mean = df_case23[df_case23['D']==1]['Y0'].mean()
    untreated_group_y0_mean = df_case23[df_case23['D']==0]['Y0'].mean()
    selection_bias = treated_group_y0_mean - untreated_group_y0_mean
    print("Explanation: The coefficient is a sum of the treatment effect and a selection bias.")
    print(f"  - The treatment effect is {tau} (Positive)")
    print(f"  - The selection bias (E[Y0|D=1] - E[Y0|D=0]) is {selection_bias:.4f} (Negative)")
    print("The large negative selection bias overwhelms the positive treatment effect, making the coefficient negative.")
    print("Conclusion for Case 2: The coefficient does NOT have to be positive.\n")


    # --- 5. Case 3: D is assigned based on X, regress Y on D and X ---
    print("--- Case 3: D is assigned based on X, regression Y ~ 1 + D + X ---")
    print("By including X in the regression, we control for the selection bias.")
    X_regr_case3 = sm.add_constant(df_case23[['D', 'X']])
    model_case3 = sm.OLS(df_case23['Y'], X_regr_case3).fit()
    coef_D_case3 = model_case3.params['D']

    print("Running regression: Y ~ 1 + D + X")
    print(f"Estimated coefficient on D: {coef_D_case3:.4f}")
    print("Explanation: By controlling for X, we eliminate the selection bias.")
    print("The coefficient on D now correctly estimates the treatment effect (which is 5.0).")
    print("Since the treatment effect is always positive, this coefficient must also be positive.")
    print("Conclusion for Case 3: The coefficient MUST be positive.\n")

    print("--- Final Summary ---")
    print("The coefficient on D must be positive in Case 1 and Case 3.")

# This ensures the script runs when executed
if __name__ == '__main__':
    run_simulation()