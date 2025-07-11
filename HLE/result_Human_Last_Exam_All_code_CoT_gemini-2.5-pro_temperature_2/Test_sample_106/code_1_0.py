import numpy as np
import statsmodels.api as sm

def demonstrate_treatment_effects():
    """
    This function simulates the three cases described in the problem
    and explains why the coefficient on the treatment D is positive or not
    based on the regression model used.
    """
    # Set a seed for reproducibility of the random data
    np.random.seed(42)

    # 1. Define a large population and their potential outcomes
    N = 100000  # Population size
    # X is the pre-treatment income
    X = np.random.uniform(20000, 80000, N)
    
    # Y0 is the potential income without treatment, which depends on X
    Y0 = 10000 + 1.2 * X + np.random.normal(0, 5000, N)

    # The individual treatment effect is positive for everyone
    individual_TE = np.random.uniform(2000, 5000, N) + 0.05 * X
    # Y1 is the potential income with treatment
    Y1 = Y0 + individual_TE
    
    print("--- Case 1: D is randomly assigned. Regress Y on D. ---")
    D1 = np.random.binomial(1, 0.5, N)
    Y_observed_1 = np.where(D1 == 1, Y1, Y0)
    
    model1 = sm.OLS(Y_observed_1, sm.add_constant(D1)).fit()
    coef_D_1 = model1.params[1]
    
    print(f"Model: Income = b0 + b1 * In_Program")
    print(f"The estimated coefficient for 'In_Program' in Case 1 is: {coef_D_1:.2f}")
    print("This coefficient is positive, as expected. With random assignment, it correctly estimates the Average Treatment Effect, which is positive by assumption.\n")

    print("--- Setup for Cases 2 & 3: D is assigned based on X (lower income individuals are more likely to be treated). ---")
    # Propensity score: P(D=1|X) decreases as X increases
    propensity = 0.9 - 0.8 * (X - X.min()) / (X.max() - X.min())
    D2 = np.random.binomial(1, propensity, N)
    Y_observed_23 = np.where(D2 == 1, Y1, Y0)
    
    print("\n--- Case 2: Conditional assignment. Regress Y on D (omitting X). ---")
    model2 = sm.OLS(Y_observed_23, sm.add_constant(D2)).fit()
    coef_D_2 = model2.params[1]

    print(f"Model: Income = b0 + b1 * In_Program")
    print(f"The estimated coefficient for 'In_Program' in Case 2 is: {coef_D_2:.2f}")
    print("This coefficient is negative. The negative selection bias (from targeting low-income individuals) is stronger than the positive treatment effect, leading to a misleading result.\n")
    
    print("--- Case 3: Conditional assignment. Regress Y on D and X. ---")
    # We create a design matrix with a constant, D2, and X
    X_for_model3 = np.column_stack((D2, X))
    model3 = sm.OLS(Y_observed_23, sm.add_constant(X_for_model3)).fit()
    coef_D_3 = model3.params[1]

    print(f"Model: Income = b0 + b1 * In_Program + b2 * Pre_Program_Income")
    print(f"The estimated coefficient for 'In_Program' in Case 3 is: {coef_D_3:.2f}")
    print("This coefficient is positive. By controlling for pre-program income (X), we remove the selection bias and recover a positive treatment effect.\n")

    print("="*60)
    print("CONCLUSION: The coefficient on D must be positive in Case 1 and Case 3.")
    print("="*60)

# Run the simulation and print the results
demonstrate_treatment_effects()

<<<E>>>