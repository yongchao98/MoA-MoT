import numpy as np
import statsmodels.api as sm

def run_simulation():
    """
    This simulation demonstrates why the coefficient on a treatment D
    is not always positive even if the true effect is, due to omitted variable bias.
    """
    # 1. Setup Simulation Parameters
    np.random.seed(42)  # for reproducibility
    N = 100_000         # number of individuals in our population
    
    # Treatment effect is a constant +5000 for everyone. It is strictly positive.
    true_treatment_effect = 5000

    print("--- Simulation Setup ---")
    print(f"Population size (N): {N}")
    print(f"True individual treatment effect (positive for everyone): {true_treatment_effect}\n")

    # 2. Generate Population Data
    # X is pre-program income. Let's assume it's between $20k and $100k.
    X = np.random.uniform(20000, 100000, N)

    # Y(0) is the potential outcome (income) without treatment.
    # It depends on pre-program income X plus some random noise.
    # Y(0) = 15000 + 1.1*X + error
    noise = np.random.normal(0, 10000, N)
    Y0 = 15000 + 1.1 * X + noise

    # Y(1) is the potential outcome with treatment.
    # Y(1) = Y(0) + Treatment Effect
    Y1 = Y0 + true_treatment_effect

    # 3. Model Treatment Assignment (The key difference between cases)
    # The jobs program (D=1) targets individuals with lower pre-program income (X).
    # We use a logistic function to model this probability.
    # P(D=1) will be high for low X and low for high X.
    logit_prob_d = 8 - 0.00015 * X
    prob_d = 1 / (1 + np.exp(-logit_prob_d))
    D = np.random.binomial(1, prob_d, N)

    # 4. Create Observed Outcome Y
    # Individuals reveal the potential outcome corresponding to their treatment status.
    Y = D * Y1 + (1 - D) * Y0

    print("--- Regression Analysis ---")
    
    # Add a constant for the intercept term in regressions
    X_with_const = sm.add_constant(X)
    D_with_const = sm.add_constant(D)
    all_vars_with_const = sm.add_constant(np.column_stack((D, X)))


    # --- Case 2 Analysis ---
    # D is conditionally random, but we regress Y on only D, omitting X.
    print("\nCase 2: Regress Y on D (omitting X)")
    model_case_2 = sm.OLS(Y, D_with_const)
    results_case_2 = model_case_2.fit()
    coef_d_case_2 = results_case_2.params[1]
    
    print("Equation: Y = {:.2f} + {:.2f}*D".format(results_case_2.params[0], coef_d_case_2))
    print("Result: The coefficient on D is NEGATIVE ({:.2f}).".format(coef_d_case_2))
    print("This is because the negative selection bias (lower-income individuals get treatment) is stronger than the positive treatment effect.")


    # --- Case 3 Analysis ---
    # D is conditionally random, and we correctly regress Y on both D and X.
    print("\nCase 3: Regress Y on D and X")
    model_case_3 = sm.OLS(Y, all_vars_with_const)
    results_case_3 = model_case_3.fit()
    coef_d_case_3 = results_case_3.params[1]

    print("Equation: Y = {:.2f} + {:.2f}*D + {:.2f}*X".format(results_case_3.params[0], coef_d_case_3, results_case_3.params[2]))
    print("Result: The coefficient on D is POSITIVE ({:.2f}).".format(coef_d_case_3))
    print("By controlling for X, the model correctly isolates the treatment effect, which is very close to the true value of 5000.")


    # --- Case 1 Hypothetical ---
    # What if D were truly random (not dependent on X)?
    print("\nHypothetical Case 1: If D were truly random")
    D_random = np.random.binomial(1, 0.5, N)
    Y_random_d = D_random * Y1 + (1-D_random) * Y0 # Recompute observed Y
    model_case_1 = sm.OLS(Y_random_d, sm.add_constant(D_random))
    results_case_1 = model_case_1.fit()
    coef_d_case_1 = results_case_1.params[1]
    
    print("Equation: Y = {:.2f} + {:.2f}*D_random".format(results_case_1.params[0], coef_d_case_1))
    print("Result: The coefficient on a truly random D is POSITIVE ({:.2f}) and also estimates the ATE.".format(coef_d_case_1))


if __name__ == '__main__':
    run_simulation()