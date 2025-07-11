import numpy as np
import statsmodels.api as sm

def run_simulations():
    """
    Simulates three scenarios to test the sign of the treatment effect coefficient.
    """
    # Set up the population/simulation
    np.random.seed(42)
    N = 100000  # Population size
    
    # X is pre-program income
    X = np.random.normal(loc=50000, scale=15000, size=N)
    
    # Potential outcome Y(0) (income without program) depends on X
    # and some random noise
    noise = np.random.normal(loc=0, scale=5000, size=N)
    Y0 = 20000 + 0.8 * X + noise

    # The treatment effect (TE) is positive for everyone
    # For simplicity, we make it a constant positive value
    TE = 2000
    
    # Potential outcome Y(1) (income with program)
    Y1 = Y0 + TE

    print("--- True Population Parameters ---")
    print(f"True Average Treatment Effect (ATE) = {np.mean(Y1 - Y0):.2f}")
    print("\n" + "="*50 + "\n")


    # --- Case 1: D is randomly assigned ---
    print("--- Case 1: D is randomly assigned. Regression: Y ~ D ---")
    # D is assigned with 50% probability
    D1 = np.random.binomial(1, 0.5, size=N)
    # Observed outcome Y
    Y_obs1 = D1 * Y1 + (1 - D1) * Y0
    # Regression model
    X_model1 = sm.add_constant(D1)
    model1 = sm.OLS(Y_obs1, X_model1).fit()
    coef_d1 = model1.params[1]
    
    print("Regression Equation: Y = B0 + B1*D")
    print(f"Estimated B0 (Constant): {model1.params[0]:.2f}")
    print(f"Estimated B1 (Coefficient on D): {coef_d1:.2f}")
    print("Conclusion: As expected, the coefficient is positive and close to the true ATE.\n")

    # --- Case 2 & 3: D is randomly assigned conditional on X ---
    print("--- Case 2 & 3: D is randomly assigned conditional on X. ---")
    # Let's assume people with lower pre-program income (X) are more likely to join
    # the program. We use a logit model for probability of treatment.
    # The term '0.0001 * (X - 50000)' makes the probability higher for X < 50000
    prob_d_cond_x = 1 / (1 + np.exp(0.0001 * (X - 50000)))
    D2 = np.random.binomial(1, prob_d_cond_x, size=N)
    # Observed outcome Y for these cases
    Y_obs2 = D2 * Y1 + (1 - D2) * Y0

    # --- Case 2: Regression Y ~ D (omitting X) ---
    print("\n--- Case 2: Regression: Y ~ D ---")
    X_model2 = sm.add_constant(D2)
    model2 = sm.OLS(Y_obs2, X_model2).fit()
    coef_d2 = model2.params[1]

    print("Regression Equation: Y = B0 + B1*D")
    print(f"Estimated B0 (Constant): {model2.params[0]:.2f}")
    print(f"Estimated B1 (Coefficient on D): {coef_d2:.2f}")
    print("Conclusion: The coefficient is negative due to strong negative selection bias.\n")

    # --- Case 3: Regression Y ~ D + X (controlling for X) ---
    print("--- Case 3: Regression: Y ~ D + X ---")
    X_model3 = sm.add_constant(np.column_stack((D2, X)))
    model3 = sm.OLS(Y_obs2, X_model3).fit()
    coef_d3 = model3.params[1]
    coef_x3 = model3.params[2]
    
    print("Regression Equation: Y = B0 + B1*D + B2*X")
    print(f"Estimated B0 (Constant): {model3.params[0]:.2f}")
    print(f"Estimated B1 (Coefficient on D): {coef_d3:.2f}")
    print(f"Estimated B2 (Coefficient on X): {coef_x3:.2f}")
    print("Conclusion: After controlling for the confounder X, the coefficient on D is positive.")

if __name__ == '__main__':
    run_simulations()
