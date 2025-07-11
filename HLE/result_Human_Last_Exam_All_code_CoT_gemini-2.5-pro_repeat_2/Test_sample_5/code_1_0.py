def calculate_proportionality_factor(d, k):
    """
    Calculates the proportionality factor C(d, k) for the expression
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C(d, k) gamma_{mu_1 ... mu_k}.

    Args:
        d (int): The number of spacetime dimensions.
        k (int): The rank of the antisymmetrized gamma matrix.

    Returns:
        int: The proportionality factor C(d, k).
    """
    if k > d:
        # The antisymmetrized product of more than d gamma matrices is zero.
        return 0
        
    # The proportionality factor is given by the formula C(d,k) = (d-k)(1-d+k)
    factor1 = d - k
    factor2 = 1 - d + k
    result = factor1 * factor2
    return result, factor1, factor2

def main():
    """
    Main function to demonstrate the calculation for example values.
    """
    # Example values for d and k
    d = 4  # Spacetime dimensions
    k = 2  # Rank of the gamma matrix product

    print(f"Calculating the proportionality factor for d = {d} and k = {k}.")
    print("The formula is C(d, k) = (d - k) * (1 - d + k).")

    factor, f1, f2 = calculate_proportionality_factor(d, k)

    print("\nThe final equation with the numbers plugged in is:")
    print(f"C({d}, {k}) = ({d} - {k}) * (1 - {d} + {k}) = ({f1}) * ({f2}) = {factor}")
    print(f"\nThe proportionality factor is: {factor}")


if __name__ == "__main__":
    main()