import sympy

def calculate_proportionality_factor(d, k):
    """
    Calculates the proportionality factor C(d, k) in the equation
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C(d, k) gamma_{mu_1 ... mu_k}
    in d-dimensional Minkowski space.

    Args:
        d (int or sympy.Symbol): The number of spacetime dimensions.
        k (int or sympy.Symbol): The rank of the antisymmetrized gamma matrix.

    Returns:
        sympy.Expr: The proportionality factor as a symbolic expression.
    """
    if not isinstance(d, sympy.Symbol):
        if not isinstance(d, int) or d < 1:
            raise ValueError("d must be a positive integer or a sympy Symbol.")
    if not isinstance(k, sympy.Symbol):
        if not isinstance(k, int) or k < 0:
            raise ValueError("k must be a non-negative integer or a sympy Symbol.")

    # The proportionality factor C(d, k)
    # C(d,k) = -((d-2k)^2 - 2k(d-k) - d)
    factor = -((d - 2*k)**2 - 2*k*(d - k) - d)
    
    return sympy.simplify(factor)

def main():
    """
    Main function to get user input and print the result.
    """
    try:
        d_str = input("Enter the number of dimensions (d): ")
        if d_str == 'd':
            d = sympy.Symbol('d')
        else:
            d = int(d_str)

        k_str = input("Enter the rank of the gamma matrix (k): ")
        if k_str == 'k':
            k = sympy.Symbol('k')
        else:
            k = int(k_str)

        factor = calculate_proportionality_factor(d, k)
        
        # Print the step-by-step thinking
        print("\nThe proportionality factor C(d,k) is given by the formula:")
        print("C(d,k) = -((d - 2k)^2 - 2k(d - k) - d)")
        
        # Print the simplified formula
        d_sym = sympy.Symbol('d')
        k_sym = sympy.Symbol('k')
        general_factor = calculate_proportionality_factor(d_sym, k_sym)
        
        print("\nSimplifying the expression, we get:")
        # Manually expand for clarity
        term1 = (d_sym - 2*k_sym)**2
        term2 = 2*k_sym*(d_sym - k_sym)
        term3 = d_sym
        expanded_inside = term1 - term2 - term3
        simplified_inside = sympy.simplify(expanded_inside)
        
        print(f"C(d,k) = -({sympy.expand(term1)} - ({sympy.expand(term2)}) - d)")
        print(f"C(d,k) = -({sympy.expand(expanded_inside)})")
        print(f"C(d,k) = {general_factor}")

        print(f"\nFor d = {d_str} and k = {k_str}, the proportionality factor is:")
        print(factor)

    except ValueError as e:
        print(f"Invalid input: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()