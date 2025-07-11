import numpy as np
import statsmodels.api as sm

def run_simulation():
    """
    Simulates three scenarios of treatment effect estimation and prints the results.
    """
    # 1. Define a population and potential outcomes
    np.random.seed(42)
    N = 100000  # Size of the population
    
    # X is pre-treatment income, normally distributed
    X = np.random.normal(loc=50000, scale=15000, size=N)
    
    # Y(0) is the income someone would have without the program.
    # It depends on their pre-treatment income X and some random noise.
    u = np.random.normal(loc=0, scale=5000, size=N)
    Y0 = 5000 + 1.0 * X + u
    
    # The treatment effect is positive for everyone. Let's say it's $2,000.
    treatment_effect = 2000
    Y1 = Y0 + treatment_effect
    
    print("--- Simulation Setup ---")
    print(f"Population size: {N}")
    print(f"True Average Treatment Effect: ${treatment_effect}")
    print("Crucial Assumption: Treatment effect is positive for everyone.")
    print("-" * 25)

    # --- Case 1: D is randomly assigned ---
    print("\n--- Case 1: D is randomly assigned. Regression: Y ~ 1 + D ---")
    D1 = np.random.binomial(1, 0.5, size=N)
    Y1_realized = Y0 * (1 - D1) + Y1 * D1
    
    # Regression
    X1_reg = sm.add_constant(D1)
    model1 = sm.OLS(Y1_realized, X1_reg).fit()
    intercept1, coef_d1 = model1.params
    
    print("The coefficient on D represents the Average Treatment Effect (ATE).")
    print("Since the true effect is positive for everyone, this coefficient must be positive.")
    print(f"Estimated Equation: Y = {intercept1:.2f} + {coef_d1:.2f}*D")
    print(f"Result: The coefficient on D is {coef_d1:.2f}, which is positive as expected.")

    # --- Case 2 & 3: D is assigned conditional on X ---
    # We model a case where the program targets lower-income individuals.
    # This will create a selection bias.
    # P(D=1|X) is higher for lower X.
    logit_p = 5 - 0.0001 * X 
    prob_d_cond_x = 1 / (1 + np.exp(-logit_p))
    D2 = np.random.binomial(1, prob_d_cond_x, size=N)
    Y23_realized = Y0 * (1 - D2) + Y1 * D2
    
    # --- Case 2: Regression Y ~ 1 + D (omitting X) ---
    print("\n--- Case 2: D depends on X. Regression: Y ~ 1 + D ---")
    X2_reg = sm.add_constant(D2)
    model2 = sm.OLS(Y23_realized, X2_reg).fit()
    intercept2, coef_d2 = model2.params
    
    print("Here, we don't control for X. The program targets lower-income people (lower X),")
    print("who would have had lower incomes (Y0) anyway. This creates a negative selection bias.")
    print("This negative bias can overwhelm the positive treatment effect.")
    print(f"Estimated Equation: Y = {intercept2:.2f} + {coef_d2:.2f}*D")
    print(f"Result: The coefficient on D is {coef_d2:.2f}, which is negative due to the strong selection bias.")

    # --- Case 3: Regression Y ~ 1 + D + X (controlling for X) ---
    print("\n--- Case 3: D depends on X. Regression: Y ~ 1 + D + X ---")
    X3_reg = sm.add_constant(np.column_stack((D2, X)))
    model3 = sm.OLS(Y23_realized, X3_reg).fit()
    intercept3, coef_d3, coef_x3 = model3.params
    
    print("By controlling for X, we eliminate the selection bias.")
    print("The coefficient on D is a weighted average of the (positive) conditional treatment effects.")
    print("Therefore, it must be positive.")
    print(f"Estimated Equation: Y = {intercept3:.2f} + {coef_d3:.2f}*D + {coef_x3:.2f}*X")
    print(f"Result: The coefficient on D is {coef_d3:.2f}, which is positive and close to the true effect.")


if __name__ == '__main__':
    run_simulation()