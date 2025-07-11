import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit

def simulate_treatment_effects():
    """
    This function simulates three scenarios to determine when the coefficient on a 
    treatment variable 'D' must be positive, given a universally positive treatment effect.
    """
    # 0. Set up a reproducible "population"
    np.random.seed(42)
    num_individuals = 100000

    # X is pre-treatment income (in $1000s)
    X = np.random.normal(loc=50, scale=10, size=num_individuals)
    
    # Random variation in income
    error = np.random.normal(loc=0, scale=5, size=num_individuals)

    # Define potential outcomes
    # Y(0) is the income without the jobs program. It depends on prior income X.
    Y0 = 5 + 1.5 * X + error
    # The treatment effect is a positive $10k for everyone.
    # Y(1) is the income with the jobs program.
    true_treatment_effect = 10
    Y1 = Y0 + true_treatment_effect

    print("This simulation demonstrates the principles using a large dataset.")
    print(f"The true treatment effect is positive for everyone: +{true_treatment_effect}\n" + "-"*60)
    
    # --- Case 1: D is randomly assigned. Regression: Y ~ D ---
    # D is assigned with a 50% probability, independent of X
    D_case1 = np.random.binomial(1, 0.5, size=num_individuals)
    Y_obs_case1 = D_case1 * Y1 + (1 - D_case1) * Y0
    
    # Run regression: Y on a constant and D
    model1 = sm.OLS(Y_obs_case1, sm.add_constant(D_case1)).fit()
    coef1 = model1.params[1]

    print("Case 1: D is randomly assigned. Regression: Y on D")
    print(f"The estimated coefficient on D is: {coef1:.4f}")
    print("Explanation: With pure random assignment, the simple regression of Y on D provides an unbiased estimate of the Average Treatment Effect (ATE). Since the treatment effect is +10 for everyone, the ATE is +10. The coefficient must be positive.")
    print("-" * 60)

    # --- Case 2: D is assigned based on X. Regression: Y ~ D ---
    # People with lower pre-treatment income X are more likely to get treatment D
    # This represents a 'selection on observables' scenario
    propensity_score = expit(5 - 0.1 * X) # `expit` is the logistic sigmoid function
    D_case2 = np.random.binomial(1, p=propensity_score, size=num_individuals)
    Y_obs_case2 = D_case2 * Y1 + (1 - D_case2) * Y0

    # Run regression: Y on a constant and D (omitting X)
    model2 = sm.OLS(Y_obs_case2, sm.add_constant(D_case2)).fit()
    coef2 = model2.params[1]
    
    print("Case 2: D is assigned based on X. Regression: Y on D")
    print(f"The estimated coefficient on D is: {coef2:.4f}")
    print("Explanation: Treatment is not random; it's correlated with X (prior income). Those with lower X are more likely to be treated. Because Y is also positively correlated with X, simply comparing the average income of treated vs. untreated is misleading. The untreated group has a higher income on average to begin with, which swamps the positive treatment effect. This creates a negative omitted variable bias. The coefficient is not guaranteed to be positive.")
    print("-" * 60)

    # --- Case 3: D is assigned based on X. Regression: Y ~ D + X ---
    # Using the same data from Case 2, but now controlling for X
    X_vars_case3 = sm.add_constant(np.stack([D_case2, X], axis=1))
    model3 = sm.OLS(Y_obs_case2, X_vars_case3).fit()
    coef3 = model3.params[1] # Coefficient for the first variable, D_case2

    print("Case 3: D is assigned based on X. Regression: Y on D and X")
    print(f"The estimated coefficient on D is: {coef3:.4f}")
    print("Explanation: By including the confounding variable X in the regression, we control for its effect. The coefficient on D now represents the treatment effect for individuals with the same level of X. Because the true effect is positive for everyone, this coefficient must also be positive and is an unbiased estimate of the true effect.")
    print("-" * 60)
    
    print("\nCONCLUSION: Based on both theory and this simulation, the coefficient on D must be positive in Case 1 and Case 3.")

simulate_treatment_effects()
<<<E>>>