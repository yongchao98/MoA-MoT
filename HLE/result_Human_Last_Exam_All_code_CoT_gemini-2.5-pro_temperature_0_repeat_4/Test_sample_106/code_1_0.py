import numpy as np
import statsmodels.api as sm

def run_simulation():
    """
    Simulates data for the three cases and runs the specified regressions
    to check the sign of the coefficient on the treatment variable D.
    """
    # Set seed for reproducibility
    np.random.seed(42)

    # 1. Data Generation Setup
    N = 100000  # Large sample size to approximate population results
    # X = pre-treatment income. Let's say it's normally distributed.
    X = np.random.normal(loc=40000, scale=10000, size=N)
    
    # Define potential outcomes
    # Y(0) = income without the program. Depends on pre-treatment income X.
    # Let's add some random noise.
    u0 = np.random.normal(0, 5000, N)
    Y0 = 20000 + 0.8 * X + u0
    
    # Y(1) = income with the program.
    # CRUCIAL ASSUMPTION: Treatment effect is positive for everyone.
    # Let's make it a constant positive effect for simplicity.
    treatment_effect = 2000
    Y1 = Y0 + treatment_effect

    print("--- Simulation Setup ---")
    print(f"Population Size (N): {N}")
    print(f"True Average Treatment Effect (ATE): {np.mean(Y1 - Y0):.2f}")
    print("Assumption: Y(1) > Y(0) for all individuals is satisfied.")
    print("-" * 25, "\n")

    # --- Case 1: D is randomly assigned ---
    print("--- Case 1: D is randomly assigned. Regress Y on D. ---")
    # D is independent of X and potential outcomes
    D1 = np.random.binomial(1, 0.5, N)
    # Observed outcome Y
    Y_obs1 = D1 * Y1 + (1 - D1) * Y0
    
    # Regression: Y ~ 1 + D
    X_sm1 = sm.add_constant(D1)
    model1 = sm.OLS(Y_obs1, X_sm1).fit()
    coef_D1 = model1.params[1]
    
    print("Regression: Y = beta_0 + beta_D * D")
    print(f"Estimated equation: Y = {model1.params[0]:.2f} + {coef_D1:.2f} * D")
    print(f"The coefficient on D is {coef_D1:.2f}.")
    print("Conclusion: As expected, the coefficient is positive, closely estimating the true ATE.")
    print("-" * 25, "\n")

    # --- Case 2 & 3: D is randomly assigned conditional on X ---
    # We will generate one dataset and use it for both Case 2 and Case 3.
    # Assignment to treatment D depends on X.
    # Let's make people with lower pre-treatment income (X) more likely to get treatment.
    # This will induce negative selection bias.
    X_centered = X - np.mean(X)
    # The probability of treatment is a logistic function of X
    # The negative sign on X_centered means lower X -> higher probability
    prob_d_equals_1 = 1 / (1 + np.exp(-(0.2 - 0.0001 * X_centered)))
    D2 = np.random.binomial(1, prob_d_equals_1, N)
    
    # Observed outcome Y for this assignment
    Y_obs2 = D2 * Y1 + (1 - D2) * Y0

    # --- Case 2: Regress Y on D (omitting X) ---
    print("--- Case 2: D is conditional on X. Regress Y on D. ---")
    # Regression: Y ~ 1 + D
    X_sm2 = sm.add_constant(D2)
    model2 = sm.OLS(Y_obs2, X_sm2).fit()
    coef_D2 = model2.params[1]
    
    print("Regression: Y = beta_0 + beta_D * D")
    print(f"Estimated equation: Y = {model2.params[0]:.2f} + ({coef_D2:.2f}) * D")
    print(f"The coefficient on D is {coef_D2:.2f}.")
    print("Conclusion: The coefficient is negative. This is because the negative selection bias")
    print("(lower-income people are selected into treatment) is stronger than the positive treatment effect.")
    print("-" * 25, "\n")

    # --- Case 3: Regress Y on D and X ---
    print("--- Case 3: D is conditional on X. Regress Y on D and X. ---")
    # Regression: Y ~ 1 + D + X
    X_sm3 = sm.add_constant(np.column_stack((D2, X)))
    model3 = sm.OLS(Y_obs2, X_sm3).fit()
    coef_D3 = model3.params[1]
    coef_X3 = model3.params[2]
    
    print("Regression: Y = beta_0 + beta_D * D + beta_X * X")
    print(f"Estimated equation: Y = {model3.params[0]:.2f} + {coef_D3:.2f} * D + {coef_X3:.2f} * X")
    print(f"The coefficient on D is {coef_D3:.2f}.")
    print("Conclusion: After controlling for the confounder X, the coefficient on D is positive,")
    print("correctly estimating the true treatment effect.")
    print("-" * 25, "\n")

if __name__ == '__main__':
    run_simulation()