import math

def calculate_proportionality_factor(d, k):
    """
    Calculates the proportionality factor for the expression
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C * gamma_{mu_1 ... mu_k}.

    Args:
        d (int): The number of spacetime dimensions.
        k (int): The rank of the antisymmetrized gamma matrix.

    Returns:
        int: The proportionality factor.
    """
    if k > d:
        print(f"The rank k={k} cannot be larger than the dimension d={d}, so the gamma matrix is zero.")
        print("The proportionality factor is trivially 0.")
        return 0
        
    # The key identity gives the factor c_k = (-1)^k * (d - 2k)
    c_k = ((-1)**k) * (d - 2 * k)
    
    # The final proportionality factor is c_k - c_k^2
    factor = c_k - c_k**2
    
    print(f"For spacetime dimension d = {d} and rank k = {k}:")
    print(f"The intermediate factor c_k is defined as (-1)^k * (d - 2k).")
    print(f"c_k = (-1)^{k} * ({d} - 2*{k}) = {c_k}")
    
    print(f"\nThe final proportionality factor is given by the formula: C = c_k - c_k^2")
    print(f"C = {c_k} - ({c_k})^2")
    print(f"C = {c_k} - {c_k**2}")
    print(f"C = {factor}")
    
    return factor

if __name__ == '__main__':
    # Example values, you can change these
    d_dimensions = 4
    k_rank = 1
    
    # Run the calculation and print the final factor
    final_factor = calculate_proportionality_factor(d_dimensions, k_rank)
    
    # For automated checking, one might uncomment the following line
    # print(f"Final Answer: {final_factor}")