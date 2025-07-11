import numpy as np
from sklearn.linear_model import LinearRegression

def run_jackknife_plus_simulation():
    """
    This function runs a simulation to demonstrate the coverage
    of the Jackknife+ prediction interval method.
    """
    # --- Step 1: Set parameters ---
    n = 50          # Number of training data points
    alpha = 0.1     # Desired miscoverage rate
    n_trials = 2000 # Number of simulation trials
    
    true_beta_0 = 2.0
    true_beta_1 = 3.0

    coverage_count = 0
    
    print(f"Running {n_trials} simulations...")
    print(f"Training set size n = {n}")
    print(f"Miscoverage rate alpha = {alpha}")
    print(f"Theoretical minimum coverage: 1 - alpha = {1 - alpha:.2f}\n")

    # --- Step 2: Main simulation loop ---
    for _ in range(n_trials):
        # --- a. Generate Data ---
        # Generate n+1 points for training data and one test point
        X_all = np.random.rand(n + 1, 1) * 10
        epsilon = np.random.randn(n + 1) * 2 # Gaussian noise
        y_all = true_beta_0 + true_beta_1 * X_all.flatten() + epsilon
        
        # The first n points are the training set D_n
        X_train, y_train = X_all[:n], y_all[:n]
        # The last point is the test point (X_{n+1}, Y_{n+1})
        X_test, y_test = X_all[n:], y_all[n:]

        # --- b. Compute LOO Scores ---
        loo_scores = []
        for i in range(n):
            # Create the leave-one-out dataset D_n \ (X_i, Y_i)
            X_loo = np.delete(X_train, i, axis=0)
            y_loo = np.delete(y_train, i)
            
            # Train algorithm A on the LOO dataset
            model_loo = LinearRegression()
            model_loo.fit(X_loo, y_loo)
            
            # Compute the LOO score on the left-out point (X_i, Y_i)
            prediction = model_loo.predict(X_train[i].reshape(1, -1))
            score = np.abs(prediction - y_train[i])[0]
            loo_scores.append(score)

        # --- c. Compute Quantile for the Prediction Interval ---
        # The set of scores S includes the n LOO scores and +infinity
        # The size of this conceptual set is n+1
        quantile_index = int(np.ceil((n + 1) * (1 - alpha)))
        
        # The quantile is the k-th smallest value. For index k, this is array element k-1
        # If quantile_index > n, the quantile would be +infinity.
        # This gives an infinitely wide interval, guaranteeing coverage.
        if quantile_index > n:
            q_hat = np.inf
        else:
            sorted_scores = sorted(loo_scores)
            q_hat = sorted_scores[quantile_index - 1]

        # --- d. Build Prediction Interval and Check Coverage ---
        # Train the full model A on the entire training set D_n
        model_full = LinearRegression()
        model_full.fit(X_train, y_train)
        
        # Get the prediction for the new test point X_{n+1}
        prediction_test = model_full.predict(X_test)[0]
        
        # Check if the true Y_{n+1} is in the prediction interval
        # The event is |prediction - y_test| <= q_hat
        if np.abs(prediction_test - y_test) <= q_hat:
            coverage_count += 1

    # --- Step 3: Report Final Results ---
    empirical_coverage = coverage_count / n_trials
    
    print("--- Simulation Results ---")
    print(f"Number of times the interval covered the true value: {coverage_count}")
    print(f"Total number of trials: {n_trials}")
    print(f"Empirical Coverage: {empirical_coverage:.4f}")
    
    lower_bound = 1-alpha
    print(f"Theoretical Lower Bound: {lower_bound:.4f}")
    if empirical_coverage >= lower_bound:
        print("Result: Empirical coverage is consistent with the theoretical guarantee (>= 1 - alpha).")
    else:
        print("Result: Empirical coverage is BELOW the theoretical guarantee.")

if __name__ == '__main__':
    run_jackknife_plus_simulation()
