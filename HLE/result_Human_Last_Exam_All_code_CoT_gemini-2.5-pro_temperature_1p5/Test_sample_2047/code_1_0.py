import numpy as np

def calculate_l(d, lambda_val):
    """
    Calculates the function l(d, λ).

    The calculation is based on the derived formula:
    l(d, λ) = (1 / (2*λ)) * [ (arccos(sqrt(2)/sqrt(d)))^2 - (arccos(sqrt(3)/sqrt(d)))^2 ]
    """
    # Check for valid inputs according to the problem statement
    if d < 4:
        print("Error: d must be greater than or equal to 4.")
        return
    if lambda_val < 1:
        print("Error: lambda must be greater than or equal to 1.")
        return

    # Calculate the arguments for arccos
    arg_for_acos_x1 = np.sqrt(3) / np.sqrt(d)
    arg_for_acos_x2 = np.sqrt(2) / np.sqrt(d)
    
    # Calculate the squared norms of the tangent vectors v(x1) and v(x2)
    norm_v1_sq = np.power(np.arccos(arg_for_acos_x1), 2)
    norm_v2_sq = np.power(np.arccos(arg_for_acos_x2), 2)

    # Calculate the final result for l(d, λ)
    result = (1 / (2 * lambda_val)) * (norm_v2_sq - norm_v1_sq)

    # Print the final equation with all the computed numbers
    print(f"For d={d} and λ={lambda_val}:")
    print(f"l(d,λ) = (1 / (2 * {lambda_val})) * ({norm_v2_sq} - {norm_v1_sq}) = {result}")

if __name__ == '__main__':
    # Set the values for d (must be >= 4) and lambda (must be >= 1)
    # The problem does not specify them, so we use the simplest valid values.
    d = 4
    lambda_val = 1
    
    calculate_l(d, lambda_val)