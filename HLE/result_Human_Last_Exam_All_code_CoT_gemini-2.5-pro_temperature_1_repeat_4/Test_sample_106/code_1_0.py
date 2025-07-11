import numpy as np
import statsmodels.api as sm

def run_simulation():
    """
    Simulates data for the three cases and prints the regression results.
    """
    # Set a seed for reproducibility
    np.random.seed(42)
    
    # 0. Define Population and Causal Model
    n_population = 100000
    # X is pre-program income, normally distributed around $40k
    X = np.random.normal(loc=40000, scale=10000, size=n_population)
    
    # Potential outcomes
    # Y(0) is the income without the program. It depends on pre-program income.
    Y0 = 20000 + 0.8 * X + np.random.normal(loc=0, scale=5000, size=n_population)
    # The treatment effect is a guaranteed positive $5,000 for everyone.
    treatment_effect = 5000
    # Y(1) is the income with the program.
    Y1 = Y0 + treatment_effect

    # 1. Case 1: D is randomly assigned
    # D is assigned to 50% of the population randomly.
    D_rand = np.random.binomial(1, 0.5, n_population)
    Y_case1 = Y0 * (1 - D_rand) + Y1 * D_rand
    
    # Run regression for Case 1: Y on a constant and D
    X_reg1 = sm.add_constant(D_rand)
    model1 = sm.OLS(Y_case1, X_reg1).fit()
    beta0_1, beta1_1 = model1.params
    
    print("--- Case 1: D is randomly assigned. Regress Y on D. ---")
    print(f"The estimated equation is: Y = {beta0_1:.2f} + {beta1_1:.2f} * D")
    print(f"The coefficient on D is {beta1_1:.2f}, which is positive as expected.\n")

    # 2. Case 2 & 3: D is assigned conditional on X
    # Let's model that people with lower pre-program income are more likely to join.
    # This creates a selection bias.
    prob_d_cond = 1 / (1 + np.exp((X - 38000) / 5000))
    D_cond = (np.random.uniform(size=n_population) < prob_d_cond).astype(int)
    Y_case23 = Y0 * (1 - D_cond) + Y1 * D_cond
    
    # Run regression for Case 2: Y on a constant and D (omitting X)
    X_reg2 = sm.add_constant(D_cond)
    model2 = sm.OLS(Y_case23, X_reg2).fit()
    beta0_2, beta1_2 = model2.params
    
    print("--- Case 2: D is conditional on X. Regress Y on D. ---")
    print(f"The estimated equation is: Y = {beta0_2:.2f} + {beta1_2:.2f} * D")
    print(f"The coefficient on D is {beta1_2:.2f}, which is negative due to omitted variable bias.\n")

    # Run regression for Case 3: Y on a constant, D, and X
    X_reg3 = sm.add_constant(np.column_stack((D_cond, X)))
    model3 = sm.OLS(Y_case23, X_reg3).fit()
    beta0_3, beta1_3, beta2_3 = model3.params
    
    print("--- Case 3: D is conditional on X. Regress Y on D and X. ---")
    print(f"The estimated equation is: Y = {beta0_3:.2f} + {beta1_3:.2f} * D + {beta2_3:.2f} * X")
    print(f"The coefficient on D is {beta1_3:.2f}, which is positive and close to the true effect.\n")

if __name__ == '__main__':
    run_simulation()