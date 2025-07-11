import numpy as np
import statsmodels.api as sm

def run_simulation():
    """
    Simulates data to demonstrate the treatment effect coefficient under three scenarios.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # 1. Define Population and Potential Outcomes
    n_samples = 100000
    # X = pre-treatment income (e.g., centered around $40k)
    X = np.random.normal(40000, 10000, n_samples)
    
    # Y(0) = potential income without treatment. Depends on pre-treatment income X.
    # We add a random error term.
    error_y0 = np.random.normal(0, 5000, n_samples)
    Y0 = 10000 + 0.8 * X + error_y0
    
    # Y(1) = potential income with treatment.
    # The treatment effect is positive for everyone. For simplicity, let's make it a constant $5,000.
    treatment_effect = 5000
    Y1 = Y0 + treatment_effect

    print("--- Simulation Setup ---")
    print(f"True Average Treatment Effect: ${treatment_effect}")
    print("Y(0) is positively correlated with X.")
    print("-" * 25 + "\n")

    # --- Case 1: D is randomly assigned ---
    print("--- Case 1: D is randomly assigned ---")
    D_case1 = np.random.binomial(1, 0.5, size=n_samples)
    Y_case1 = Y0 * (1 - D_case1) + Y1 * D_case1
    X1_model = sm.add_constant(D_case1)
    model1 = sm.OLS(Y_case1, X1_model).fit()
    print("Regression: Y ~ 1 + D")
    print(f"Estimated coefficient on D: {model1.params[1]:.2f}")
    print("Result: As expected, the coefficient is positive and close to the true ATE.\n")

    # --- Case 2 & 3: D is assigned based on X ---
    # People with lower pre-treatment income X are more likely to get treatment D.
    # This will create a negative selection bias.
    logit_p = (np.mean(X) - X) / 5000  # Higher probability for lower X
    prob_d = 1 / (1 + np.exp(-logit_p))
    D_case23 = np.random.binomial(1, prob_d, size=n_samples)
    Y_case23 = Y0 * (1 - D_case23) + Y1 * D_case23
    
    # Check for selection bias
    avg_y0_treated = np.mean(Y0[D_case23 == 1])
    avg_y0_untreated = np.mean(Y0[D_case23 == 0])
    selection_bias = avg_y0_treated - avg_y0_untreated

    # --- Case 2: Omitted variable bias ---
    print("--- Case 2: D depends on X, X is omitted from regression ---")
    print(f"Selection Bias (E[Y(0)|D=1] - E[Y(0)|D=0]): {selection_bias:.2f}")
    X2_model = sm.add_constant(D_case23)
    model2 = sm.OLS(Y_case23, X2_model).fit()
    print("Regression: Y ~ 1 + D")
    print(f"Estimated coefficient on D: {model2.params[1]:.2f}")
    print("Result: The coefficient is negative because the negative selection bias is stronger than the positive treatment effect.\n")

    # --- Case 3: Controlling for the confounder ---
    print("--- Case 3: D depends on X, X is included in regression ---")
    X3_model = sm.add_constant(np.column_stack((D_case23, X)))
    model3 = sm.OLS(Y_case23, X3_model).fit()
    print("Regression: Y ~ 1 + D + X")
    print(f"Estimated coefficient on D: {model3.params[1]:.2f}")
    print("Result: After controlling for X, the coefficient on D is positive and close to the true ATE.")

run_simulation()