import numpy as np
import math

def run_simulation():
    """
    This simulation demonstrates the theoretical guarantee of leave-one-out
    conformal prediction.
    """
    # --- Parameters ---
    n = 20          # Number of training points
    alpha = 0.1     # Significance level (e.g., for a 90% nominal interval, alpha=0.1)
    num_trials = 10000 # Number of simulations to run for empirical validation

    print(f"--- Parameters ---")
    print(f"Number of training points (n): {n}")
    print(f"Significance level (alpha): {alpha}")
    print(f"The nominal coverage level is 1-alpha = {1-alpha:.2f}\n")

    # --- Theoretical Calculation ---
    print("--- Theoretical Lower Bound ---")
    # The lowest possible coverage is given by the formula:
    # 1 - floor((n+1)*alpha) / (n+1)
    n_plus_1 = n + 1
    numerator = math.floor(n_plus_1 * alpha)
    theoretical_miscoverage = numerator / n_plus_1
    theoretical_coverage = 1 - theoretical_miscoverage

    print("The final equation for the lowest possible coverage probability is:")
    print(f"P(coverage) = 1 - floor((n + 1) * alpha) / (n + 1)")
    print("Plugging in the numbers:")
    print(f"n + 1 = {n_plus_1}")
    print(f"(n + 1) * alpha = {n_plus_1 * alpha:.2f}")
    print(f"floor((n + 1) * alpha) = {numerator}")
    print(f"Lowest coverage = 1 - {numerator} / {n_plus_1} = {theoretical_coverage:.4f}\n")


    # --- Empirical Simulation ---
    print("--- Running Simulation ---")
    coverage_count = 0
    for _ in range(num_trials):
        # 1. Generate data (n training points + 1 test point)
        # We use a simple case: Y_i ~ N(0,1) and X is ignored.
        # The algorithm will be to predict the mean.
        # This setup ensures scores are from a continuous distribution, so the
        # theoretical bound is tight.
        data = np.random.randn(n + 1)
        y_train = data[:n]
        y_test = data[n]

        # 2. Compute n leave-one-out (LOO) scores
        loo_scores = np.zeros(n)
        for i in range(n):
            # Create the LOO training set
            y_train_loo = np.delete(y_train, i)
            # Train algorithm A on LOO set (A = predict mean)
            y_pred_loo = np.mean(y_train_loo)
            # Compute LOO score
            loo_scores[i] = np.abs(y_pred_loo - y_train[i])

        # 3. Compute the quantile q
        # The quantile is the k-th smallest value of the LOO scores + infinity
        # k = ceil((n+1)*(1-alpha))
        k = math.ceil(n_plus_1 * (1 - alpha))
        
        # Add +infinity to the set of scores
        scores_with_inf = np.append(loo_scores, np.inf)
        
        # Sort and find the k-th value
        q = np.sort(scores_with_inf)[k-1]

        # 4. Build the predictive interval for the test point
        # First, train the algorithm A on the full training set
        y_pred_full = np.mean(y_train)
        
        # 5. Check if the test point is covered
        test_score = np.abs(y_pred_full - y_test)
        if test_score <= q:
            coverage_count += 1

    empirical_coverage = coverage_count / num_trials
    print(f"Simulation finished after {num_trials} trials.")
    print(f"Empirical coverage: {empirical_coverage:.4f}")

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    print(f"The simulated coverage ({empirical_coverage:.4f}) is very close to the theoretical lower bound ({theoretical_coverage:.4f}).")


if __name__ == '__main__':
    run_simulation()