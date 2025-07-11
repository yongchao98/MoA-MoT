from fractions import Fraction

def get_bernoulli_numbers():
    """
    Returns a dictionary of the first few even-indexed Bernoulli numbers.
    Convention B_1 = -1/2. The B_{2k} are standard.
    """
    return {
        2: Fraction(1, 6),
        4: Fraction(-1, 30),
        6: Fraction(1, 42)
    }

def zeta_neg_int(n):
    """
    Computes the Riemann zeta function at a negative integer -n, for odd n.
    The formula is zeta(1-k) = -B_k/k. So zeta(-n) = zeta(1-(n+1)) = -B_{n+1}/(n+1).
    n must be odd.
    """
    k = n + 1
    B = get_bernoulli_numbers()
    if k not in B:
        raise ValueError(f"Bernoulli number B_{k} is not available in the list.")
    return -B[k] / k

def solve():
    """
    Calculates the orbifold Euler characteristic of the quotient stack [U/G].
    """
    # Genus of a plane quartic curve is g=3.
    g = 3

    print("Step 1: Identify the space and the formula.")
    print("The stack of smooth plane quartics modulo PGL(3) is the moduli stack of non-hyperelliptic genus 3 curves, M_3^nh.")
    print("The orbifold Euler characteristic is calculated as chi(M_3^nh) = chi(M_3) - chi(H_3).\n")

    # Calculate chi(M_3)
    print("Step 2: Calculate chi(M_3) using the Harer-Zagier formula chi(M_g) = zeta(1-2g).")
    B = get_bernoulli_numbers()
    k_m3 = 2 * g
    zeta_arg_m3 = 1 - k_m3
    
    chi_M3 = -B[k_m3] / k_m3
    print(f"For g=3, chi(M_3) = zeta({zeta_arg_m3}) = -B_{k_m3}/{k_m3}.")
    print(f"B_{k_m3} = {B[k_m3]}")
    print(f"chi(M_3) = -({B[k_m3]}) / {k_m3} = {chi_M3}\n")

    # Calculate chi(H_3)
    print("Step 3: Calculate chi(H_3) using the Bergström-Tommasi formula.")
    print("For g=3, the formula is chi(H_3) = (1/2)*zeta(-5) - (1/2)*zeta(-1)*zeta(-3).\n")
    
    # Values needed for H_3
    zeta_m1 = zeta_neg_int(1)
    k_m1 = 2
    zeta_m3 = zeta_neg_int(3)
    k_m3_zeta = 4
    zeta_m5 = chi_M3 # zeta(-5)
    
    print("Calculating required zeta function values:")
    print(f"zeta(-1) = -B_{k_m1}/{k_m1} = -({B[k_m1]})/{k_m1} = {zeta_m1}")
    print(f"zeta(-3) = -B_{k_m3_zeta}/{k_m3_zeta} = -({B[k_m3_zeta]})/{k_m3_zeta} = {zeta_m3}")
    print(f"zeta(-5) = {zeta_m5}\n")

    term1 = Fraction(1, 2) * zeta_m5
    term2 = Fraction(1, 2) * zeta_m1 * zeta_m3
    chi_H3 = term1 - term2
    
    print(f"chi(H_3) = (1/2) * ({zeta_m5}) - (1/2) * ({zeta_m1}) * ({zeta_m3})")
    print(f"         = {term1} - ({term2})")
    print(f"         = {chi_H3}\n")
    
    # Final calculation
    result = chi_M3 - chi_H3
    print("Step 4: Calculate the final result.")
    print(f"chi(M_3^nh) = chi(M_3) - chi(H_3)")
    print(f"            = {chi_M3} - ({chi_H3})")
    print(f"            = {result}\n")

    print(f"The final answer is {result.numerator}/{result.denominator}.")

solve()
<<< -47/20160 >>>