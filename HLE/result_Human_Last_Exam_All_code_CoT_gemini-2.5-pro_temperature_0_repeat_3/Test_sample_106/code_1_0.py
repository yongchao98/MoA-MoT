import numpy as np
import statsmodels.api as sm
from scipy.special import expit

def run_simulation():
    """
    Simulates data for the three cases and runs the specified regressions.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # 1. Define a population
    n_population = 100000
    # Generate pre-treatment income X
    X = np.random.uniform(20000, 90000, n_population)

    # 2. Define potential outcomes Y(0) and Y(1)
    # Y(0) is the income without the program, depends on pre-treatment income X
    # Add some random noise
    error_y0 = np.random.normal(0, 5000, n_population)
    Y0 = 10000 + 1.1 * X + error_y0

    # Y(1) is the income with the program.
    # The treatment effect (TE) is Y(1) - Y(0). We ensure it's always positive.
    treatment_effect = 2000 + 0.05 * X + np.random.uniform(500, 1000, n_population)
    Y1 = Y0 + treatment_effect

    # --- Case 1: D is randomly assigned ---
    print("--- Case 1: D is randomly assigned. Regress Y on D. ---")
    # D is a coin flip, independent of X
    D1 = np.random.binomial(1, 0.5, n_population)
    # Observed outcome Y
    Y_case1 = D1 * Y1 + (1 - D1) * Y0
    # Regression: Y on a constant and D
    X_case1 = sm.add_constant(D1)
    model1 = sm.OLS(Y_case1, X_case1).fit()
    print(f"The equation is: Y = {model1.params[0]:.2f} + {model1.params[1]:.2f}*D")
    print(f"The coefficient on D is {model1.params[1]:.2f}, which is positive.\n")


    # --- Case 2 & 3: D assignment depends on X ---
    # Let's model that people with lower pre-treatment income X are more likely to join the program
    # We use a logistic function for the probability
    prob_d_equals_1 = expit(8 - 0.00015 * X)
    D2 = np.random.binomial(1, prob_d_equals_1, n_population)
    # Observed outcome Y for these cases
    Y_case23 = D2 * Y1 + (1 - D2) * Y0

    # --- Case 2: Regress Y on D (omitting X) ---
    print("--- Case 2: D assignment depends on X. Regress Y on D. ---")
    X_case2 = sm.add_constant(D2)
    model2 = sm.OLS(Y_case23, X_case2).fit()
    print(f"The equation is: Y = {model2.params[0]:.2f} + {model2.params[1]:.2f}*D")
    print(f"The coefficient on D is {model2.params[1]:.2f}. Due to selection bias, it is negative in this simulation, so it does not have to be positive.\n")

    # --- Case 3: Regress Y on D and X ---
    print("--- Case 3: D assignment depends on X. Regress Y on D and X. ---")
    # Create the design matrix with a constant, D, and X
    X_case3 = sm.add_constant(np.column_stack((D2, X)))
    model3 = sm.OLS(Y_case23, X_case3).fit()
    print(f"The equation is: Y = {model3.params[0]:.2f} + {model3.params[1]:.2f}*D + {model3.params[2]:.2f}*X")
    print(f"The coefficient on D is {model3.params[1]:.2f}, which is positive as the regression now controls for X.\n")

if __name__ == '__main__':
    run_simulation()