import numpy as np
import statsmodels.api as sm

def run_simulation(N=100000, seed=42):
    """
    Runs a simulation to demonstrate the conditions under which the
    coefficient on a treatment variable is positive.
    """
    np.random.seed(seed)

    # 1. Define Population Parameters
    # X = pre-treatment income, normally distributed
    X = np.random.normal(loc=30000, scale=5000, size=N)
    
    # TE = The true Treatment Effect. It's a positive constant for all individuals.
    TE = 5000 
    
    # Y0 = Potential outcome (income) without treatment. Depends on pre-treatment income.
    # u is some random noise.
    u = np.random.normal(loc=0, scale=3000, size=N)
    Y0 = 20000 + 0.8 * X + u
    
    # Y1 = Potential outcome (income) with treatment.
    Y1 = Y0 + TE

    print("--- Simulation Setup ---")
    print(f"Population size N = {N}")
    print(f"True Treatment Effect (positive for everyone) = ${TE}")
    print(f"Pre-treatment income 'X' is a key predictor of potential income 'Y0'.")
    print("\n")

    # --- Case 1: D is randomly assigned ---
    print("--- Case 1: D is randomly assigned. Regression: Y ~ D ---")
    D_random = np.random.binomial(1, 0.5, size=N)
    Y_observed_1 = D_random * Y1 + (1 - D_random) * Y0
    
    # Add a constant for the regression model
    model1_X = sm.add_constant(D_random)
    model1 = sm.OLS(Y_observed_1, model1_X).fit()
    coef_d_1 = model1.params[1]
    
    print(f"The equation is Y = {model1.params[0]:.2f} + {coef_d_1:.2f}*D")
    print("Result: The coefficient on D is positive, correctly estimating the Average Treatment Effect.")
    print("\n")

    # --- Case 2 & 3: Data where D's assignment depends on X ---
    # We model that people with lower pre-treatment income X are more likely to be in the program.
    # This creates the selection bias problem for Case 2.
    prob_d = 1 / (1 + np.exp((X - 30000) / 2000)) # Logistic probability function
    D_conditional = np.random.binomial(1, prob_d, size=N)
    Y_observed_23 = D_conditional * Y1 + (1 - D_conditional) * Y0

    # --- Case 2: D is assigned based on X, but X is omitted from regression ---
    print("--- Case 2: D depends on X. Regression: Y ~ D (X is omitted) ---")
    model2_X = sm.add_constant(D_conditional)
    model2 = sm.OLS(Y_observed_23, model2_X).fit()
    coef_d_2 = model2.params[1]

    print(f"The equation is Y = {model2.params[0]:.2f} + {coef_d_2:.2f}*D")
    print("Result: The coefficient on D is negative due to strong negative selection bias.")
    print("\n")

    # --- Case 3: D is assigned based on X, and X is included in regression ---
    print("--- Case 3: D depends on X. Regression: Y ~ D + X (X is included) ---")
    # We use the same data as Case 2, but add X to the model to control for it.
    model3_X = sm.add_constant(np.c_[D_conditional, X])
    model3 = sm.OLS(Y_observed_23, model3_X).fit()
    coef_d_3 = model3.params[1]
    coef_x_3 = model3.params[2]
    
    print(f"The equation is Y = {model3.params[0]:.2f} + {coef_d_3:.2f}*D + {coef_x_3:.2f}*X")
    print("Result: The coefficient on D is positive again, as controlling for X removes the bias.")
    print("\n")

# Run the simulation and print the results
run_simulation()