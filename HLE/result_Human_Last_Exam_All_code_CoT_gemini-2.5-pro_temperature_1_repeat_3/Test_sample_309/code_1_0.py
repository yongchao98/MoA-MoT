def solve_bundle_adjustment_problem(N, M):
    """
    Calculates the maximum number of landmarks that can be marginalized in BA.

    Args:
        N (int): The number of landmarks.
        M (int): The number of cameras.

    Returns:
        int: The maximum number of landmarks that can be marginalized.
    """
    # To make the bundle adjustment problem solvable with a unique solution,
    # we need to fix the 7-DOF gauge freedom (3 rotation, 3 translation, 1 scale).
    # A minimal way to define a 3D coordinate system is by using 3 landmarks.
    # These 3 landmarks become the fixed reference frame and are not part of the optimization variables.
    # The remaining N - 3 landmarks are variables and can be marginalized using the Schur complement.
    
    if N < 3:
        print("Warning: At least 3 landmarks are required to define a stable reference frame.")
        # In this case, the problem is ill-defined, but for the formula, we return 0.
        max_marginalized_landmarks = 0
        equation_str = f"Result = 0 (since N < 3)"
    else:
        max_marginalized_landmarks = N - 3
        # We need to output each number in the final equation.
        # The equation is: N - 3 = result
        num1 = N
        num2 = 3
        result = max_marginalized_landmarks
        equation_str = f"Maximum landmarks to marginalize = {num1} - {num2} = {result}"
        
    print(f"Given N = {N} landmarks and M = {M} cameras.")
    print(equation_str)
    
    return max_marginalized_landmarks

# Example values, you can change them
N = 100
M = 10

# Calculate and print the result
solve_bundle_adjustment_problem(N, M)
