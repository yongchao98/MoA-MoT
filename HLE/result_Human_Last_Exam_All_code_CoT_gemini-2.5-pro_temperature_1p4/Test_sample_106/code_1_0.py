import numpy as np
import statsmodels.api as sm

def run_simulation():
    """
    Simulates data for a job training program to analyze treatment effect coefficients.

    We will create a population where a job training program (D) always has a
    positive effect on future income (Y), but people with lower prior income (X)
    are more likely to enter the program.
    """
    # Set parameters for the simulation
    np.random.seed(42)
    num_individuals = 100000
    true_treatment_effect = 2.0  # The program increases everyone's income by $2,000

    # 1. Generate population data
    # X = pre-program income (in $1000s), normally distributed
    X = np.random.normal(loc=30, scale=5, size=num_individuals)
    # Unobserved factors affecting income
    noise = np.random.normal(loc=0, scale=2, size=num_individuals)

    # 2. Define potential outcomes
    # Y(0) is income without the program. It depends on prior income X.
    Y0 = 5 + 1.0 * X + noise
    # Y(1) is income with the program. It's always higher than Y(0) for everyone.
    Y1 = Y0 + true_treatment_effect

    # --- Case 1: D is randomly assigned ---
    # Each person has a 50% chance of being in the program, regardless of X
    D1 = np.random.binomial(1, 0.5, size=num_individuals)
    Y_obs1 = Y1 * D1 + Y0 * (1 - D1)
    
    # Regression for Case 1: Y on D
    X_model1 = sm.add_constant(D1)
    model1 = sm.OLS(Y_obs1, X_model1).fit()
    coef_case1 = model1.params[1]

    # --- Case 2 & 3: D is assigned based on X ---
    # People with lower pre-program income (X) are more likely to be treated.
    # We use a logistic function to model the probability of treatment.
    # The probability of D=1 decreases as X increases.
    logit_prob = 2 - 0.1 * X 
    prob_d = 1 / (1 + np.exp(-logit_prob))
    D2 = np.random.binomial(1, prob_d, size=num_individuals)
    Y_obs2 = Y1 * D2 + Y0 * (1 - D2)

    # --- Regression for Case 2: Y on D (omitting X) ---
    X_model2 = sm.add_constant(D2)
    model2 = sm.OLS(Y_obs2, X_model2).fit()
    coef_case2 = model2.params[1]
    
    # --- Regression for Case 3: Y on D and X ---
    X_model3 = sm.add_constant(np.column_stack((D2, X)))
    model3 = sm.OLS(Y_obs2, X_model3).fit()
    coef_case3 = model3.params[1]

    # Print the results
    print("This simulation demonstrates the three cases.")
    print(f"The true treatment effect is positive and constant: {true_treatment_effect}\n")

    print(f"--- Case 1: D is randomly assigned. Regress Y on D. ---")
    print(f"The estimated coefficient on D is: {coef_case1:.4f}")
    print("This is very close to the true treatment effect, as expected.\n")
    
    print(f"--- Case 2: D depends on X. Regress Y on D. ---")
    print(f"The estimated coefficient on D is: {coef_case2:.4f}")
    print("The coefficient is negative! This is due to omitted variable bias.")
    print("People with low pre-program income (X) are selected into the program,")
    print("and also have lower post-program income (Y) on average, creating a negative bias.\n")

    print(f"--- Case 3: D depends on X. Regress Y on D and X. ---")
    print(f"The estimated coefficient on D is: {coef_case3:.4f}")
    print("By controlling for X, we remove the bias and recover the positive treatment effect.\n")

if __name__ == '__main__':
    run_simulation()
<<<E>>>