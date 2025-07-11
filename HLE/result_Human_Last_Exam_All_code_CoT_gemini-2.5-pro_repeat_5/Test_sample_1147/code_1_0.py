import numpy as np
from sklearn.linear_model import LinearRegression
import math

def simulate_loo_coverage(n, alpha, num_trials=2000):
    """
    Simulates the Leave-One-Out conformal prediction process to verify coverage.

    Args:
        n (int): Number of training points.
        alpha (float): Significance level (e.g., 0.1 for 90% confidence).
        num_trials (int): Number of simulation trials to run.
    """
    coverage_count = 0

    for _ in range(num_trials):
        # 1. Generate data (n+1 points)
        # We use a simple linear model with noise. The method works for any model/data.
        X = np.random.rand(n + 1, 1) * 10
        noise = np.random.randn(n + 1) * 2
        y = 2 * X.squeeze() + 5 + noise

        X_train, y_train = X[:n], y[:n]
        X_test, y_test = X[n:], y[n:]

        # 2. Compute LOO scores
        loo_scores = np.zeros(n)
        for i in range(n):
            # Create the leave-one-out dataset
            X_loo = np.delete(X_train, i, axis=0)
            y_loo = np.delete(y_train, i)

            # Train algorithm A on the LOO set
            model_loo = LinearRegression()
            model_loo.fit(X_loo, y_loo)

            # Compute the score for the left-out point
            pred_i = model_loo.predict(X_train[i:i+1])
            loo_scores[i] = np.abs(pred_i - y_train[i])

        # 3. Build the predictive interval
        # Train final model on all n points
        model_final = LinearRegression()
        model_final.fit(X_train, y_train)

        # Calculate the quantile Q
        # The set of scores is S = {s_1, ..., s_n, +inf}
        # The quantile q_{1-alpha}(S) is the k-th smallest element where k = ceil((1-alpha)*(n+1))
        # This corresponds to np.quantile with interpolation='higher'
        scores_with_inf = np.append(loo_scores, np.inf)
        q = np.quantile(scores_with_inf, 1 - alpha, method='higher')

        # 4. Check coverage for the test point
        y_pred_test = model_final.predict(X_test)
        test_score = np.abs(y_pred_test - y_test)

        if test_score <= q:
            coverage_count += 1

    # Calculate probabilities
    empirical_coverage = coverage_count / num_trials
    
    # The exact theoretical probability is ceil((1-alpha)*(n+1)) / (n+1)
    k = math.ceil((1 - alpha) * (n + 1))
    theoretical_coverage = k / (n + 1)
    
    lower_bound = 1 - alpha

    print(f"--- Simulation for n={n}, alpha={alpha} ---")
    print(f"Empirical coverage over {num_trials} trials: {empirical_coverage:.4f}")
    print(f"Exact theoretical coverage: {theoretical_coverage:.4f}")
    print(f"Guaranteed lower bound (1 - alpha): {lower_bound:.4f}")
    print("-" * 35 + "\n")

# --- Case 1: (1-alpha)*(n+1) is an integer ---
# Here, the theoretical coverage should be exactly 1-alpha.
n1 = 19
alpha1 = 0.05  # (1-0.05)*(19+1) = 0.95 * 20 = 19 (integer)
simulate_loo_coverage(n1, alpha1)

# --- Case 2: (1-alpha)*(n+1) is not an integer ---
# Here, the theoretical coverage should be strictly greater than 1-alpha.
n2 = 20
alpha2 = 0.05 # (1-0.05)*(20+1) = 0.95 * 21 = 19.95 (not integer)
simulate_loo_coverage(n2, alpha2)

# The question asks for the lowest value P can take, which is 1-alpha.
# The code should output the numbers in this final equation.
final_alpha = 0.1 # Using a generic alpha for the final print
final_result = 1 - final_alpha
print("The final equation for the lowest possible coverage probability is:")
print(f"1 - alpha = result")
print(f"1 - {final_alpha} = {final_result}")
