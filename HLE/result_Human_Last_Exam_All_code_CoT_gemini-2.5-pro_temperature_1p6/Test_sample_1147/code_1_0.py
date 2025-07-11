import numpy as np
from sklearn.linear_model import LinearRegression

def demonstrate_loo_coverage(n, alpha, num_simulations=2000):
    """
    Runs a simulation to numerically verify the coverage property of
    the leave-one-out conformal prediction method.

    Args:
        n (int): Number of training points.
        alpha (float): Desired significance level.
        num_simulations (int): Number of Monte Carlo trials to run.
    """
    coverage_count = 0
    # Data generating process: Y = 2*X + 5 + noise
    # Algorithm: Linear Regression

    for _ in range(num_simulations):
        # 1. Generate data D_n and a test point (X_{n+1}, Y_{n+1})
        X_all = np.random.rand(n + 1, 1) * 10
        epsilon = np.random.randn(n + 1) # Standard normal noise
        Y_all = (2 * X_all.flatten() + 5 + epsilon).flatten()

        # Split into training data D_n and test point
        X_train, Y_train = X_all[:n], Y_all[:n]
        X_test, Y_test = X_all[n:], Y_all[n:]

        # 2. Compute the n leave-one-out (LOO) scores
        loo_scores = np.zeros(n)
        for i in range(n):
            # Create LOO dataset D_n \ (X_i, Y_i)
            X_loo = np.delete(X_train, i, axis=0)
            Y_loo = np.delete(Y_train, i, axis=0)

            # Train algorithm on the LOO set
            model_loo = LinearRegression()
            model_loo.fit(X_loo, Y_loo)

            # Compute the LOO score (residual) for the held-out point (X_i, Y_i)
            prediction_i = model_loo.predict(X_train[i].reshape(1, -1))
            loo_scores[i] = np.abs(prediction_i - Y_train[i])

        # 3. Train the final model A on the full training set D_n
        model_full = LinearRegression()
        model_full.fit(X_train, Y_train)

        # 4. Build the prediction interval for the test point X_{n+1}
        # The set of scores includes the n LOO scores and +infinity
        S = np.append(loo_scores, [np.inf])
        
        # The quantile level k is determined by alpha
        k = int(np.ceil((n + 1) * (1 - alpha)))
        
        # The interval width is the k-th smallest score
        q_value = np.sort(S)[k - 1]

        # Make prediction for the new point
        prediction_test = model_full.predict(X_test)
        
        # 5. Check for coverage: is Y_{n+1} inside the interval?
        # Y_test is in [pred - q, pred + q] <=> |pred - Y_test| <= q
        test_residual = np.abs(prediction_test - Y_test)
        if test_residual <= q_value:
            coverage_count += 1
            
    empirical_coverage = coverage_count / num_simulations
    # The exact theoretical coverage for distributions with no ties
    theoretical_coverage = np.ceil((n + 1) * (1 - alpha)) / (n + 1)

    print("--- Leave-One-Out Coverage Simulation ---")
    print(f"Parameters: n = {n}, alpha = {alpha}")
    print(f"Empirical coverage from {num_simulations} simulations: {empirical_coverage:.4f}")
    print(f"Theoretical coverage guarantee: P(coverage) >= {1-alpha}")
    print(f"Exact theoretical coverage for this case: {theoretical_coverage:.4f}\n")
    
# --- Main Answer ---
# The analysis shows the coverage is guaranteed to be at least 1-alpha.
# This bound is tight (i.e., achievable).
# We now print the final answer in the requested format.
alpha_val = 0.1 # A chosen alpha for the demonstration
# To make the theoretical coverage exactly 1-alpha, we can choose n=9.
# Then (n+1)*(1-alpha) = 10*0.9 = 9, which is an integer.
demonstrate_loo_coverage(n=9, alpha=alpha_val)

print("--- Final Answer ---")
print("The lowest value that P(Y_{n+1} in C_{n+1}(X_{n+1})) can take is described by the equation:")
print("Lowest Value = 1 - \u03B1") # u03B1 is the unicode for alpha
print(f"\nFor any given alpha, this defines the greatest lower bound on the coverage probability.")
print(f"For example, if alpha were {alpha_val}, the numbers in the final equation would be:")
print(f"Lowest Value = 1 - {alpha_val} = {1 - alpha_val}")
