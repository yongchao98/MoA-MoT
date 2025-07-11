import numpy as np
import statsmodels.api as sm

def run_simulation():
    """
    Simulates three scenarios to test the sign of a treatment effect coefficient.
    """
    # 1. Setup a population
    np.random.seed(42)
    n_population = 100000
    # X = pre-program income, from $20k to $100k
    X = np.random.uniform(20000, 100000, n_population)
    # Potential outcome without treatment (Y(0)) depends on pre-program income X
    # plus some random noise
    noise = np.random.normal(0, 5000, n_population)
    Y0 = 15000 + 1.1 * X + noise
    
    # The treatment effect is ALWAYS positive for everyone
    treatment_effect = 5000
    Y1 = Y0 + treatment_effect

    print("--- Simulation Setup ---")
    print(f"Population size: {n_population}")
    print(f"True Treatment Effect for everyone: +{treatment_effect}")
    print("\n--- Running Regressions ---\n")

    # Case 1: D is randomly assigned
    # D is a coin flip, independent of X and Y
    D_case1 = np.random.randint(0, 2, n_population)
    Y_case1 = Y0 * (1 - D_case1) + Y1 * D_case1
    X1 = sm.add_constant(D_case1)
    model1 = sm.OLS(Y_case1, X1).fit()
    print("Case 1: D is random, Regress Y on D")
    print(f"Estimated coefficient on D: {model1.params[1]:.2f}")
    print("Result: The coefficient is positive, as expected.\n")

    # Case 2 & 3: D is randomly assigned conditional on X
    # Probability of treatment is higher for lower pre-program income X
    # This creates a negative correlation between D and X
    prob_d = 1 / (1 + np.exp(-12 + 0.0002 * X))
    D_case23 = (np.random.uniform(0, 1, n_population) < prob_d).astype(int)
    Y_case23 = Y0 * (1 - D_case23) + Y1 * D_case23

    # Case 2: Regress Y on D (omitting X)
    X2 = sm.add_constant(D_case23)
    model2 = sm.OLS(Y_case23, X2).fit()
    print("Case 2: D depends on X, Regress Y on D (omitting X)")
    print(f"Estimated coefficient on D: {model2.params[1]:.2f}")
    print("Result: The coefficient is negative. The negative selection bias is stronger than the positive treatment effect.\n")

    # Case 3: Regress Y on D and X
    X3 = sm.add_constant(np.column_stack((D_case23, X)))
    model3 = sm.OLS(Y_case23, X3).fit()
    print("Case 3: D depends on X, Regress Y on D and X")
    # The coefficient is at index 1 because index 0 is the constant and index 2 is X
    print(f"Estimated coefficient on D: {model3.params[1]:.2f}")
    print(f"Estimated coefficient on X: {model3.params[2]:.2f}")
    print("Result: The coefficient on D is positive, successfully recovering the true effect after controlling for X.")

run_simulation()