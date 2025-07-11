import numpy as np
from sklearn.linear_model import LinearRegression

def run_conformal_simulation(n, alpha, n_trials=5000):
    """
    Simulates the leave-one-out conformal prediction process to verify coverage probability.

    Args:
        n (int): Number of training points.
        alpha (float): Miscoverage level (e.g., 0.1 for 90% confidence).
        n_trials (int): Number of Monte Carlo simulations to run.
    """
    coverage_count = 0
    
    # Store the last interval for demonstration
    last_interval_info = {}

    print(f"--- Running Simulation ---")
    print(f"n = {n}, alpha = {alpha}, trials = {n_trials}")
    print(f"Target coverage level (1 - alpha): {1 - alpha:.4f}")
    
    # Theoretical coverage P = ceil((n+1)(1-alpha)) / (n+1)
    k_numerator = np.ceil((n + 1) * (1 - alpha))
    theoretical_coverage = k_numerator / (n + 1)
    print(f"Theoretical exact coverage: {k_numerator}/{n+1} = {theoretical_coverage:.4f}")

    for i in range(n_trials):
        # 1. Generate data (n training points + 1 test point)
        # Y = 2*X + 5 + noise
        X_full = np.random.rand(n + 1, 1) * 10
        noise = np.random.randn(n + 1, 1) * 2
        Y_full = 2 * X_full + 5 + noise
        
        X_train, Y_train = X_full[:n], Y_full[:n]
        X_test, Y_test = X_full[n:], Y_full[n:]

        # 2. Compute LOO scores
        loo_scores = []
        for j in range(n):
            # Create the LOO training set by removing point j
            X_loo = np.delete(X_train, j, axis=0)
            Y_loo = np.delete(Y_train, j, axis=0)
            
            # Train algorithm A (Linear Regression) on the LOO set
            model_loo = LinearRegression()
            model_loo.fit(X_loo, Y_loo)
            
            # Compute the score for the left-out point
            y_pred_loo = model_loo.predict(X_train[j:j+1])
            score = np.abs(y_pred_loo - Y_train[j])
            loo_scores.append(score[0,0])

        # 3. Get the quantile q
        # Add +infinity to the set of scores
        scores_with_inf = np.append(loo_scores, np.inf)
        
        # Calculate the quantile index k
        # We use k_numerator which was already calculated
        # Note: k_numerator is 1-based, so we subtract 1 for 0-based index
        quantile_idx = int(k_numerator - 1)
        
        # The quantile q is the k-th smallest score
        q = np.sort(scores_with_inf)[quantile_idx]

        # 4. Train final model on all n training points
        model_final = LinearRegression()
        model_final.fit(X_train, Y_train)
        
        # 5. Build the predictive interval for the test point
        y_pred_final = model_final.predict(X_test)[0, 0]
        interval_min = y_pred_final - q
        interval_max = y_pred_final + q

        # 6. Check for coverage
        if interval_min <= Y_test[0, 0] <= interval_max:
            coverage_count += 1
            
        # Save last run's info
        if i == n_trials - 1:
            last_interval_info = {
                "y_hat": y_pred_final,
                "q": q,
                "interval_min": interval_min,
                "interval_max": interval_max
            }

    # Calculate and print final results
    empirical_coverage = coverage_count / n_trials
    print(f"\n--- Results ---")
    print(f"Empirical Coverage: {empirical_coverage:.4f}")
    
    # Print the details of the final equation from the last trial
    info = last_interval_info
    print("\n--- Example Prediction Interval (from last trial) ---")
    print(f"C(X_test) = [y_hat - q, y_hat + q]")
    # Using the print function to show each number in the equation
    print(f"Equation values: y_hat = {info['y_hat']:.4f}, q = {info['q']:.4f}")
    print(f"Final Interval: [{info['interval_min']:.4f}, {info['interval_max']:.4f}]")


# Case 1: (n+1)(1-alpha) is an integer. Coverage should be exactly 1-alpha.
run_conformal_simulation(n=19, alpha=0.1)

print("\n" + "="*40 + "\n")

# Case 2: (n+1)(1-alpha) is not an integer. Coverage should be > 1-alpha.
run_conformal_simulation(n=20, alpha=0.1)
