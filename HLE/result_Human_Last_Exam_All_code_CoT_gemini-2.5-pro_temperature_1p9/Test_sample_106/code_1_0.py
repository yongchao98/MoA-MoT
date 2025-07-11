import numpy as np
import statsmodels.api as sm

def run_treatment_effect_simulation():
    """
    Simulates three scenarios of treatment effect regression to determine
    when the coefficient on the treatment variable must be positive.
    """
    # Step 1: Set up a simulated population
    np.random.seed(42)  # for reproducibility
    num_individuals = 200000

    # X is pre-program income, drawn from a log-normal distribution
    X = np.random.lognormal(mean=10.5, sigma=0.5, size=num_individuals)

    # Step 2: Define potential outcomes
    # Y(0) is the income without the jobs program. It depends on pre-program income.
    Y0 = 1.0 * X + np.random.normal(0, 500, num_individuals)
    # The treatment effect is a constant positive amount, ensuring Y(1) > Y(0) for everyone.
    # The problem states the effect is positive for everyone. We set it to $3,000.
    treatment_effect = 3000
    Y1 = Y0 + treatment_effect

    print("--- Simulation Setup ---")
    print(f"The true treatment effect is a positive constant: ${treatment_effect}\n")

    # --- Case 1: D is randomly assigned ---
    # Treatment is assigned with a 50% probability, independent of X.
    D_case1 = np.random.binomial(1, 0.5, size=num_individuals)
    # Observed outcome Y is determined by treatment status D.
    Y_case1 = D_case1 * Y1 + (1 - D_case1) * Y0

    # Regress Y on D
    model1_X = sm.add_constant(D_case1)
    model1 = sm.OLS(Y_case1, model1_X).fit()
    coef_D_case1 = model1.params[1]

    print("--- Case 1: D is randomly assigned. Regression: Y ~ D ---")
    print("In this case, D is independent of all other factors. The coefficient on D should estimate the Average Treatment Effect.")
    print(f"The estimated coefficient on D is: {coef_D_case1:.2f}")
    if coef_D_case1 > 0:
        print("Result: The coefficient is positive, as expected.\n")
    else:
        print("Result: The coefficient is NOT positive.\n")


    # --- Case 2 & 3: D is randomly assigned conditional on X ---
    # The probability of receiving treatment depends on pre-program income X.
    # Individuals with lower income are more likely to be in the program.
    # This introduces a selection bias.
    prob_d = 1 / (1 + np.exp((X - np.median(X)) / 10000))
    D_case23 = np.random.binomial(1, prob_d, size=num_individuals)
    Y_case23 = D_case23 * Y1 + (1 - D_case23) * Y0

    # Demonstrate the selection bias
    avg_X_treated = X[D_case23 == 1].mean()
    avg_X_untreated = X[D_case23 == 0].mean()
    print("--- Selection Bias Check for Cases 2 & 3 ---")
    print(f"Average pre-income (X) for the treated group: ${avg_X_treated:.2f}")
    print(f"Average pre-income (X) for the untreated group: ${avg_X_untreated:.2f}")
    print("The treated group has a much lower average pre-income, creating a strong selection bias.\n")

    # --- Case 2: Conditional assignment, short regression ---
    # Regress Y on only D, omitting the confounding variable X.
    model2_X = sm.add_constant(D_case23)
    model2 = sm.OLS(Y_case23, model2_X).fit()
    coef_D_case2 = model2.params[1]
    
    print("--- Case 2: D assigned based on X. Regression: Y ~ D ---")
    print("Here we omit X from the regression. The coefficient on D will be biased.")
    print("The selection bias (lower pre-income people get treated) can overwhelm the positive treatment effect.")
    print(f"The estimated coefficient on D is: {coef_D_case2:.2f}")
    if coef_D_case2 > 0:
        print("Result: The coefficient is positive.\n")
    else:
        print("Result: The coefficient is NEGATIVE, demonstrating this case is not guaranteed to be positive.\n")

    # --- Case 3: Conditional assignment, long regression ---
    # Regress Y on both D and X, controlling for the confounder.
    model3_X = sm.add_constant(np.c_[D_case23, X])
    model3 = sm.OLS(Y_case23, model3_X).fit()
    coef_D_case3 = model3.params[1]

    print("--- Case 3: D assigned based on X. Regression: Y ~ D + X ---")
    print("By including X in the regression, we control for the selection bias.")
    print("The coefficient on D should now estimate the true positive treatment effect.")
    print(f"The estimated coefficient on D is: {coef_D_case3:.2f}")
    if coef_D_case3 > 0:
        print("Result: The coefficient is positive, as expected.\n")
    else:
        print("Result: The coefficient is NOT positive.\n")

    print("--- Conclusion ---")
    print("The coefficient on D was positive in Case 1 and Case 3, but negative in Case 2.")
    print("Therefore, only in Case 1 and Case 3 must the coefficient on D be positive in the population.")

if __name__ == '__main__':
    run_treatment_effect_simulation()