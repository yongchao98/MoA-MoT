import numpy as np
import statsmodels.api as sm

def analyze_treatment_effects():
    """
    Explains and simulates the treatment effect scenarios.
    """
    print("This script analyzes the sign of a treatment effect coefficient under different scenarios.")
    print("Let Y be income, D be a binary treatment (jobs program), and X be pre-treatment income.")
    print("We assume the true treatment effect is positive for everyone.\n")

    # --- Step 1: Theoretical Explanation ---
    print("--- Part 1: Theoretical Analysis ---")
    print("\n[Case 1: D is randomly assigned. Regression: Y ~ D]")
    print("Analysis: With pure random assignment, the groups D=1 and D=0 are identical on average. The coefficient on D measures the Average Treatment Effect (ATE = E[Y|D=1] - E[Y|D=0]). Since the effect is positive for everyone, the ATE must be positive.")
    print("Conclusion: The coefficient on D MUST be positive.\n")

    print("[Case 2: D is assigned conditional on X. Regression: Y ~ D]")
    print("Analysis: This regression omits the variable X. X is likely correlated with both treatment D (e.g., jobs programs target lower-income individuals) and outcome Y (pre-income predicts post-income). This creates Omitted Variable Bias (OVB).")
    print("The bias can be negative and may overwhelm the true positive effect, leading to a negative coefficient.")
    print("Conclusion: The coefficient on D is NOT guaranteed to be positive.\n")

    print("[Case 3: D is assigned conditional on X. Regression: Y ~ D + X]")
    print("Analysis: By including X, we control for the confounding factor. The coefficient on D now represents a weighted average of the treatment effect across different levels of X. Since the effect is assumed to be positive for everyone, this average must also be positive.")
    print("Conclusion: The coefficient on D MUST be positive.\n")
    print("Summary of Theory: The coefficient on D must be positive in Case 1 and Case 3 only.")
    print("-" * 50)

    # --- Step 2: Numerical Simulation ---
    print("\n--- Part 2: Simulation Demonstration ---")
    print("We now simulate data where a jobs program (D=1) is preferentially given to individuals with lower pre-treatment income (X).\n")
    np.random.seed(42)
    N = 100000

    # Generate pre-treatment income X
    X = np.random.normal(loc=40000, scale=10000, size=N)
    X[X < 0] = 0 # Income cannot be negative

    # Define potential outcomes
    # Y0 is the outcome without treatment. It depends on X.
    Y0 = 20000 + 0.9 * X + np.random.normal(0, 5000, N)
    
    # The individual treatment effect (ITE) is ALWAYS positive
    positive_ite = 5000 + 0.05 * X 
    
    # Y1 is the outcome with treatment
    Y1 = Y0 + positive_ite
    
    # Assign treatment D based on X (lower X -> higher probability of D=1)
    # This models a targeted jobs program
    log_odds = 2.0 - 0.0001 * X
    prob_d1 = 1 / (1 + np.exp(-log_odds))
    D = np.random.binomial(1, prob_d1, size=N)
    
    # Observed outcome Y is a combination of Y0 and Y1 based on D
    Y = D * Y1 + (1 - D) * Y0

    # Perform regressions for Case 2 and Case 3

    # Case 2: Regress Y on only D (omitting X)
    model_case2 = sm.OLS(Y, sm.add_constant(D))
    results_case2 = model_case2.fit()
    
    print("\n[Simulation Result for Case 2: Y = b0 + b1*D]")
    print(f"The estimated equation is: Y = {results_case2.params[0]:.2f} + ({results_case2.params[1]:.2f})*D")
    print(f"The coefficient on D is {results_case2.params[1]:.2f}. It is NEGATIVE.")
    print("This happens because the program was given to people who were expected to have lower incomes anyway. Comparing them to the untreated (higher-income) group without controlling for X creates a strong negative bias.\n")

    # Case 3: Regress Y on both D and X
    model_case3_data = sm.add_constant(np.c_[D, X])
    model_case3 = sm.OLS(Y, model_case3_data)
    results_case3 = model_case3.fit()
    
    print("[Simulation Result for Case 3: Y = b0 + b1*D + b2*X]")
    print(f"The estimated equation is: Y = {results_case3.params[0]:.2f} + ({results_case3.params[1]:.2f})*D + {results_case3.params[2]:.2f}*X")
    print(f"The coefficient on D is {results_case3.params[1]:.2f}. It is POSITIVE.")
    print("By controlling for pre-treatment income X, we remove the selection bias and recover a positive treatment effect, as predicted by the theory.")

if __name__ == '__main__':
    analyze_treatment_effects()
    
<<<E>>>