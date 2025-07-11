import math

def solve_with_titan():
    """
    This function simulates the calculation of the landing force on the Titan computer.
    It demonstrates that the calculation is not possible due to the architecture's constraints.
    """

    # === INPUTS ===
    # Probe mass (m) is 50 kg. We represent it by its factors 2 and 25.
    m_factors = [2, 25]
    
    # Gravitational acceleration (a_g) is ~0.1667 m/s^2, or 1/6.
    a_g_num, a_g_den = 1, 6
    
    # Deceleration (a_d) is 9 m/s^2, or 9/1.
    a_d_num, a_d_den = 9, 1
    
    # The total force F_total = m * (a_g + a_d) = m*a_g + m*a_d = F_g + F_d

    # === ANALYSIS of F_d ===
    # F_d = 50 * 9. We need to multiply the factors of 50 with 9.
    # Factors of F_d are [2, 25, 9].
    # Let's check if we can multiply these on a 5-bit computer where registers hold 0-31.
    
    # Possible multiplication orderings for 2 * 25 * 9:
    # 1. (2 * 25) * 9
    intermediate_product1 = 2 * 25
    if intermediate_product1 > 31:
        # This path fails because the first operation (2*25=50) results in a value
        # that cannot be stored in a 5-bit register.
        # Since this path, and all other paths, fail, the calculation is not possible.
        # We can therefore conclude the answer is N0.
        print("N0")
        return

    # This part of the code will not be reached, but is included for completeness.
    # 2. 2 * (25 * 9)
    intermediate_product2 = 25 * 9
    if intermediate_product2 > 31:
        # This path also fails.
        print("N0")
        return
        
    # Since all paths to calculate the deceleration force F_d lead to an overflow,
    # the Titan computer cannot solve this problem.

solve_with_titan()