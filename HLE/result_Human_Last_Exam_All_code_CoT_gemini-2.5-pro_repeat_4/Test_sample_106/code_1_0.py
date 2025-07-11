import numpy as np
import statsmodels.api as sm

def run_simulation(N=100000):
    """
    Simulates data and runs regressions for the three cases.
    """
    print("Running simulation to find the sign of the coefficient on D in each case...")
    print("-" * 70)

    # --- Data Generating Process (DGP) Setup ---
    # We define a DGP where X (pre-treatment income) affects both Y (post-treatment income)
    # and the probability of treatment D. This will create omitted variable bias in Case 2.
    
    # Let's assume a strong positive relationship between pre- and post-income
    # and a constant, positive treatment effect for everyone.
    true_ate = 5.0
    beta_x_strong = 3.0
    
    # Generate X (pre-treatment income) from a uniform distribution
    np.random.seed(42)
    X = np.random.uniform(low=10, high=40, size=N)
    
    # Define potential outcomes
    # Y(0) is the income someone would have without treatment.
    # Y(1) is the income someone would have with treatment.
    # The treatment effect is Y(1) - Y(0) = true_ate = 5.0, which is positive for everyone.
    noise = np.random.normal(0, 5, size=N)
    Y_0 = 20 + beta_x_strong * X + noise
    Y_1 = Y_0 + true_ate

    # --- Case 1: D is randomly assigned ---
    # D is independent of X. Let's assume a 50% chance of treatment.
    D_case1 = np.random.binomial(1, 0.5, size=N)
    Y_case1 = D_case1 * Y_1 + (1 - D_case1) * Y_0
    
    # Regression: Y on a constant and D
    X_reg_case1 = sm.add_constant(D_case1)
    model_case1 = sm.OLS(Y_case1, X_reg_case1).fit()
    coef_D_case1 = model_case1.params[1]
    
    print("Case 1: D is randomly assigned. Regress Y on D.")
    print(f"The estimated equation is: Y = {model_case1.params[0]:.2f} + {coef_D_case1:.2f} * D")
    print(f"The coefficient on D is {coef_D_case1:.4f}. It is positive as expected, close to the true ATE of {true_ate}.")
    print("-" * 70)

    # --- Case 2 & 3: D is randomly assigned conditional on X ---
    # We make treatment (D) more likely for individuals with lower pre-treatment income (X).
    # This creates a negative correlation between D and X.
    prob_d_equals_1 = 1 / (1 + np.exp(0.2 * X - 5)) # A logistic function
    D_case23 = np.random.binomial(1, prob_d_equals_1, size=N)
    Y_case23 = D_case23 * Y_1 + (1 - D_case23) * Y_0
    
    # --- Case 2: Regress Y on D (omitting X) ---
    X_reg_case2 = sm.add_constant(D_case23)
    model_case2 = sm.OLS(Y_case23, X_reg_case2).fit()
    coef_D_case2 = model_case2.params[1]
    
    print("Case 2: D is assigned based on X. Regress Y on D.")
    print(f"The estimated equation is: Y = {model_case2.params[0]:.2f} + {coef_D_case2:.2f} * D")
    print(f"The coefficient on D is {coef_D_case2:.4f}. It is negative due to omitted variable bias.")
    print("-" * 70)
    
    # --- Case 3: Regress Y on D and X (controlling for X) ---
    X_reg_case3 = sm.add_constant(np.column_stack((D_case23, X)))
    model_case3 = sm.OLS(Y_case23, X_reg_case3).fit()
    coef_D_case3 = model_case3.params[1]

    print("Case 3: D is assigned based on X. Regress Y on D and X.")
    print(f"The estimated equation is: Y = {model_case3.params[0]:.2f} + {coef_D_case3:.2f} * D + {model_case3.params[2]:.2f} * X")
    print(f"The coefficient on D is {coef_D_case3:.4f}. It is positive, recovering the true effect.")
    print("-" * 70)
    
    # Final conclusion from simulation
    print("Summary:")
    print(f"Case 1 (Y ~ D, random assignment): Coefficient on D is {'positive' if coef_D_case1 > 0 else 'not positive'}.")
    print(f"Case 2 (Y ~ D, assignment on X): Coefficient on D is {'positive' if coef_D_case2 > 0 else 'not positive'}.")
    print(f"Case 3 (Y ~ D + X, assignment on X): Coefficient on D is {'positive' if coef_D_case3 > 0 else 'not positive'}.")

# Run the analysis
run_simulation()