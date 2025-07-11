import math

def solve_hss_submatrices():
    """
    Calculates the number of submatrices accessed in an HSS matrix
    multiplication for a given tree depth.
    """
    # The depth of the Hierarchical Semi-separable tree
    depth = 4

    # The formula to calculate the number of submatrices is 3 * 2^(d-1) - 2.
    # This is derived from:
    # - Number of leaf nodes (diagonal blocks): 2^(d-1)
    # - Number of non-leaf nodes: Sum of 2^k for k=0 to d-2, which is 2^(d-1) - 1
    # - Number of off-diagonal blocks: 2 * (number of non-leaf nodes)
    # - Total = (2 * (2^(d-1) - 1)) + 2^(d-1) = 3 * 2^(d-1) - 2

    # Perform the calculation step-by-step
    power_val = depth - 1
    power_result = int(math.pow(2, power_val))
    mult_result = 3 * power_result
    final_result = mult_result - 2

    # Print the explanation and the final equation step-by-step
    print(f"To find the number of submatrices accessed in an HSS tree of depth d = {depth}, we use the formula: 3 * 2^(d-1) - 2.")
    print("\nHere is the step-by-step calculation:")
    print(f"Equation: 3 * 2^({depth}-1) - 2")
    print(f"Step 1: Calculate the exponent -> 3 * 2^{power_val} - 2")
    print(f"Step 2: Calculate the power of 2 -> 3 * {power_result} - 2")
    print(f"Step 3: Perform the multiplication -> {mult_result} - 2")
    print(f"Step 4: Perform the subtraction -> {final_result}")
    print(f"\nSo, for a depth of {depth}, the final equation is: 3 * {power_result} - 2 = {final_result}")

solve_hss_submatrices()