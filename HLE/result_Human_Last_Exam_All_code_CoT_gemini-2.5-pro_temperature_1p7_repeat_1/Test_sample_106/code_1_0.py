import numpy as np
import statsmodels.api as sm
import pandas as pd

def run_simulation():
    """
    Simulates data to test the three cases and prints the regression results.
    """
    # Set seed for reproducibility
    np.random.seed(42)

    # 1. DATA GENERATION
    # Let's create a population of N individuals
    N = 200000

    # X is pre-treatment income. To make the confounding clear, let's create two groups:
    # half with low pre-treatment income (X=10k), half with high (X=50k).
    X = np.random.choice([10, 50], N, p=[0.5, 0.5])

    # Y(0) is the income without the jobs program. It depends on pre-treatment income.
    # We add some random noise.
    error = np.random.normal(0, 5, N)
    Y0 = 15 + 0.8 * X + error

    # Y(1) is the income with the jobs program. The problem states the treatment effect is positive for everyone.
    # Let's say the program adds a constant $5k to everyone's income.
    treatment_effect = 5
    Y1 = Y0 + treatment_effect

    print("--- Simulation Setup ---")
    print(f"True Treatment Effect: {treatment_effect}")
    print("This effect is positive for all individuals by construction (Y1 > Y0).\n")


    # --- CASE 1: D is randomly assigned ---
    print("--- CASE 1: Y ~ D (D is Unconditionally Random) ---")
    # Treatment is assigned to 50% of the population, independent of X.
    D_case1 = np.random.binomial(1, 0.5, N)
    Y_case1 = D_case1 * Y1 + (1 - D_case1) * Y0
    
    # Regression: Y on a constant and D
    X_reg_case1 = sm.add_constant(D_case1)
    model_case1 = sm.OLS(Y_case1, X_reg_case1).fit()
    coef_d_case1 = model_case1.params[1]

    print(f"Regression model: Y = b0 + b1*D")
    print(f"Estimated equation: Y = {model_case1.params[0]:.2f} + {coef_d_case1:.2f}*D")
    print(f"Conclusion: The coefficient on D is {coef_d_case1:.2f}, which is positive as expected.")
    print("This estimates the Average Treatment Effect (ATE).\n")


    # --- CASES 2 & 3: D is randomly assigned conditional on X ---
    # People with lower pre-treatment income (X) are more likely to be in the program.
    # P(D=1 | X=10) = 0.8; P(D=1 | X=50) = 0.2
    prob_d = np.where(X == 10, 0.8, 0.2)
    D_case23 = np.random.binomial(1, prob_d, N)
    Y_case23 = D_case23 * Y1 + (1 - D_case23) * Y0

    # --- CASE 2: Short regression, omitting X ---
    print("--- CASE 2: Y ~ D (D is Conditionally Random, X omitted) ---")
    # Regression: Y on a constant and D
    X_reg_case2 = sm.add_constant(D_case23)
    model_case2 = sm.OLS(Y_case23, X_reg_case2).fit()
    coef_d_case2 = model_case2.params[1]
    
    print(f"Regression model: Y = b0 + b1*D")
    print(f"Estimated equation: Y = {model_case2.params[0]:.2f} + {coef_d_case2:.2f}*D")
    print(f"Conclusion: The coefficient on D is {coef_d_case2:.2f}, which is NEGATIVE.")
    print("This is due to omitted variable bias. The untreated group has a much higher pre-treatment income on average, leading to a higher average outcome, which masks the positive treatment effect.\n")

    # --- CASE 3: Long regression, controlling for X ---
    print("--- CASE 3: Y ~ D + X (D is Conditionally Random, X included) ---")
    # Regression: Y on a constant, D, and X
    df = pd.DataFrame({'Y': Y_case23, 'D': D_case23, 'X': X})
    X_reg_case3 = sm.add_constant(df[['D', 'X']])
    model_case3 = sm.OLS(df['Y'], X_reg_case3).fit()
    coef_d_case3 = model_case3.params['D']

    print(f"Regression model: Y = b0 + b1*D + b2*X")
    print(f"Estimated equation: Y = {model_case3.params['const']:.2f} + {coef_d_case3:.2f}*D + {model_case3.params['X']:.2f}*X")
    print(f"Conclusion: The coefficient on D is {coef_d_case3:.2f}, which is positive as expected.")
    print("By controlling for X, we remove the confounding and correctly estimate the positive treatment effect.\n")


# Run the simulation and print results
run_simulation()