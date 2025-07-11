import numpy as np
import statsmodels.api as sm

def run_causal_simulation():
    """
    Simulates three scenarios to test when the treatment effect coefficient is positive.
    Y = income, D = in jobs program, X = pre-treatment income.
    We assume the true treatment effect is positive for everyone.
    """
    # 1. Setup: Define a population and potential outcomes
    N = 200000  # Number of individuals in the population
    np.random.seed(42)  # For reproducibility

    # Pre-treatment income X, from $20k to $80k
    X = np.random.uniform(20000, 80000, N)

    # Potential outcomes: Y(0) is income without program, Y(1) is income with program.
    # We model Y as a function of pre-treatment income X.
    Y0 = 15000 + 1.1 * X
    # The treatment effect is a constant +$5000 for everyone. This ensures Y(1) > Y(0).
    true_treatment_effect = 5000
    Y1 = Y0 + true_treatment_effect

    print("--- Simulation Setup ---")
    print(f"The true treatment effect is positive for everyone: ${true_treatment_effect}")
    print("\n" + "="*60 + "\n")

    # 2. Case 1: D is randomly assigned (like a pure RCT)
    print("--- Case 1: D is randomly assigned ---")
    print("Model: Y ~ 1 + D")
    D_case1 = np.random.binomial(1, 0.5, size=N)  # 50% get the program, irrespective of X
    Y_obs_case1 = D_case1 * Y1 + (1 - D_case1) * Y0
    
    # Run the regression
    model1 = sm.OLS(Y_obs_case1, sm.add_constant(D_case1)).fit()
    coef_D1 = model1.params[1]
    
    print(f"The estimated coefficient on D is: {coef_D1:.2f}")
    print("Result: The coefficient is positive, correctly estimating the ATE.")
    print("Reason: With random assignment, there is no confounding. E[Y|D=1]-E[Y|D=0] = ATE.")
    print("\n" + "="*60 + "\n")

    # 3. Setup for Cases 2 & 3: D is assigned based on X
    # Let's assume the program targets lower-income individuals.
    # Probability of treatment decreases as pre-treatment income X increases.
    prob_D_conditional = 1.25 - (X / 80000)
    D_case23 = np.random.binomial(1, prob_D_conditional)
    Y_obs_case23 = D_case23 * Y1 + (1 - D_case23) * Y0

    # This creates confounding: treated people have lower X on average.
    mean_X_treated = np.mean(X[D_case23 == 1])
    mean_X_untreated = np.mean(X[D_case23 == 0])
    print(f"--- Setup for Cases 2 & 3: D depends on X ---")
    print(f"Average pre-treatment income (X) for Treated (D=1): ${mean_X_treated:,.2f}")
    print(f"Average pre-treatment income (X) for Untreated (D=0): ${mean_X_untreated:,.2f}")
    print("This difference in X creates confounding (omitted variable bias).")
    print("\n" + "="*60 + "\n")

    # 4. Case 2: Conditional assignment, but regression omits X
    print("--- Case 2: D is assigned conditional on X, but X is omitted from regression ---")
    print("Model: Y ~ 1 + D")
    
    model2 = sm.OLS(Y_obs_case23, sm.add_constant(D_case23)).fit()
    coef_D2 = model2.params[1]

    print(f"The estimated coefficient on D is: {coef_D2:.2f}")
    print("Result: The coefficient is negative!")
    print("Reason: Omitted variable bias. The negative effect of having a lower starting income (X) is stronger than the positive treatment effect.")
    print("\n" + "="*60 + "\n")

    # 5. Case 3: Conditional assignment, and regression includes X
    print("--- Case 3: D is assigned conditional on X, and X is included in regression ---")
    print("Model: Y ~ 1 + D + X")
    
    X_vars_case3 = np.column_stack((D_case23, X))
    model3 = sm.OLS(Y_obs_case23, sm.add_constant(X_vars_case3)).fit()
    coef_D3 = model3.params[1]

    print(f"The estimated coefficient on D is: {coef_D3:.2f}")
    print("Result: The coefficient is positive, correctly estimating the treatment effect.")
    print("Reason: By controlling for X, we have removed the confounding bias.")
    print("\n" + "="*60 + "\n")

    print("Final Conclusion: The coefficient on D must be positive in Case 1 and Case 3.")

run_causal_simulation()
<<<E>>>