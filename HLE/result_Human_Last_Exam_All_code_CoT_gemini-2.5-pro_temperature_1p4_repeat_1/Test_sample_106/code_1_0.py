import numpy as np
import statsmodels.api as sm

def demonstrate_treatment_effects():
    """
    This function simulates data for the three cases to demonstrate
    when the coefficient on a treatment variable must be positive.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # 1. Define the population and Data Generating Process (DGP)
    N = 100_000  # Population size
    # X is pre-treatment income, normally distributed around a mean of 40 (in thousands)
    X = np.random.normal(loc=40, scale=10, size=N)
    
    # Potential outcome without treatment (Y0) depends on pre-treatment income X
    # plus some random individual variation.
    Y0 = 20 + 1.1 * X + np.random.normal(loc=0, scale=5, size=N)
    
    # The treatment effect (TE) is positive for everyone, a random value between 2 and 5.
    TE = np.random.uniform(low=2, high=5, size=N)
    
    # Potential outcome with treatment (Y1) is Y0 + TE.
    # This guarantees Y1 > Y0 for all individuals.
    Y1 = Y0 + TE
    
    print("--- Simulation Setup ---")
    print(f"Average True Treatment Effect: {np.mean(TE):.4f}\n")

    # --- CASE 1: D is randomly assigned. Regress Y on D. ---
    print("--- Case 1: D is randomly assigned ---")
    # Treatment is assigned to 50% of the population randomly.
    D1 = np.random.binomial(1, 0.5, size=N)
    # Observed outcome Y depends on whether the individual received the treatment.
    Y_obs1 = D1 * Y1 + (1 - D1) * Y0
    
    # Run regression of Y on a constant and D.
    X1_reg = sm.add_constant(D1)
    model1 = sm.OLS(Y_obs1, X1_reg).fit()
    coef_d1 = model1.params[1]
    
    print("Regression: Y ~ 1 + D")
    print(f"The coefficient on D is: {coef_d1:.4f}")
    print("Result: The coefficient is positive, correctly estimating the Average Treatment Effect (ATE).\n")


    # --- CASE 2 & 3: D is randomly assigned conditional on X ---
    # Assignment probability depends on X. People with lower income (X) are more likely to be treated.
    # This creates selection bias if X is not controlled for.
    logit_p = 4.0 - 0.1 * X
    prob_d = 1 / (1 + np.exp(-logit_p))
    D23 = np.random.binomial(1, prob_d, size=N)
    Y_obs23 = D23 * Y1 + (1 - D23) * Y0

    # --- CASE 2: Regress Y on D (omitting X) ---
    print("--- Case 2: D depends on X, but we only regress Y ~ 1 + D ---")
    X2_reg = sm.add_constant(D23)
    model2 = sm.OLS(Y_obs23, X2_reg).fit()
    coef_d2 = model2.params[1]

    print("Regression: Y ~ 1 + D")
    print(f"The coefficient on D is: {coef_d2:.4f}")
    print("Result: The coefficient is NEGATIVE. This is due to selection bias.")
    print("The control group has higher pre-income (X) and thus higher Y on average, hiding the positive treatment effect.\n")

    # --- CASE 3: Regress Y on D and X ---
    print("--- Case 3: D depends on X, and we regress Y ~ 1 + D + X ---")
    # We use the same data as Case 2, but now include X in the regression.
    X3_reg = sm.add_constant(np.column_stack((D23, X)))
    model3 = sm.OLS(Y_obs23, X3_reg).fit()
    coef_d3 = model3.params[1] # Coefficient for D is the second one

    print("Regression: Y ~ 1 + D + X")
    print(f"The coefficient on D is: {coef_d3:.4f}")
    print("Result: The coefficient is positive. By controlling for X, we remove the selection bias and recover the treatment effect.")

# Run the simulation and print the results
demonstrate_treatment_effects()
<<<E>>>