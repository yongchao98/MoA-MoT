import math

def calculate_imitation_learning_performance_bound(H, A, lambda_val, R_max=1.0):
    """
    Calculates the tightest upper bound on the performance difference 
    J(pi^*) - J(pi_hat) based on a standard imitation learning result.

    Args:
        H (int): The episode horizon.
        A (int): The size of the action space, |A|.
        lambda_val (float): A hyperparameter of the algorithm.
        R_max (float, optional): The maximum possible reward. Defaults to 1.0.

    Returns:
        float: The calculated upper bound for J(pi^*) - J(pi_hat).
    """
    # The population total variation (TV) risk is bounded by: |A| * (1 - exp(-lambda))
    tv_risk_bound = A * (1 - math.exp(-lambda_val))
    
    # The performance difference is bounded by: (H * (H + 1) / 2) * R_max * TV_risk
    # This is a standard result from the analysis of compounding errors in imitation learning.
    compounding_error_coeff = (H * (H + 1)) / 2
    
    performance_bound = compounding_error_coeff * R_max * tv_risk_bound
    
    return performance_bound

# The problem is symbolic, so no specific values for H, |A|, or lambda are provided.
# We will use example values to demonstrate the calculation and show the numbers in the equation.
H_example = 20
A_example = 5
lambda_example = 0.5
R_max_example = 1.0 # A common assumption when not specified.

# Calculate the bound for the example values
final_bound = calculate_imitation_learning_performance_bound(H_example, A_example, lambda_example, R_max_example)

# Print the general formula and the result for the example case, showing each number.
print("The tightest upper bound for J(pi^*) - J(pi_hat) is given by the formula:")
print("Bound = (H * (H + 1) / 2) * R_max * |A| * (1 - exp(-lambda))")
print("\n--- Example Calculation ---")
print(f"Using example parameters: H = {H_example}, |A| = {A_example}, lambda = {lambda_example}, R_max = {R_max_example}")
print(f"1. The compounding error coefficient is (H * (H + 1) / 2):")
print(f"   ({H_example} * ({H_example} + 1) / 2) = {(H_example * (H_example + 1)) / 2}")
print(f"2. The upper bound on the TV risk is |A| * (1 - exp(-lambda)):")
print(f"   {A_example} * (1 - exp(-{lambda_example})) = {A_example * (1 - math.exp(-lambda_example))}")
print(f"3. The final upper bound on the performance difference is:")
print(f"   Result = {((H_example * (H_example + 1)) / 2)} * {R_max_example} * {A_example * (1 - math.exp(-lambda_example))}")
print(f"   J(pi^*) - J(pi_hat) <= {final_bound}")
