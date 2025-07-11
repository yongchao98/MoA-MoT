import math
from sympy import bernoulli
from fractions import Fraction

def solve_euler_characteristic():
    """
    Calculates the orbifold Euler characteristic of the moduli space of smooth plane quartics.
    """
    # Genus of a plane quartic curve d=4 is g = (d-1)(d-2)/2
    g = 3
    
    print("The problem asks for the orbifold Euler characteristic of the quotient stack [U/G],")
    print("where U is the space of smooth plane quartic curves and G = PGL(3).")
    print("This stack is the moduli stack of non-hyperelliptic genus 3 curves, M_3^nh.")
    print("\nWe use the decomposition of the moduli stack of genus 3 curves, M_3, into its")
    print("non-hyperelliptic (M_3^nh) and hyperelliptic (M_3^h) loci:")
    print("chi_orb(M_3^nh) = chi_orb(M_3) - chi_orb(M_3^h)\n")

    # Step 1: Calculate chi_orb(M_3)
    print("--- Step 1: Compute chi_orb(M_3) ---")
    print("The Harer-Zagier formula for the orbifold Euler characteristic of M_g is:")
    print(f"chi_orb(M_g) = -B_2g / (4*g)")
    
    # We need the 2g-th Bernoulli number, B_6
    b_2g_idx = 2 * g
    # sympy.bernoulli returns a sympy.Rational
    B_2g_sympy = bernoulli(b_2g_idx)
    B_2g = Fraction(B_2g_sympy.p, B_2g_sympy.q)
    
    print(f"For g = {g}, we need the Bernoulli number B_{b_2g_idx} = {B_2g}.")
    
    chi_orb_M3_num = -B_2g.numerator
    chi_orb_M3_den = B_2g.denominator * 4 * g
    chi_orb_M3 = Fraction(chi_orb_M3_num, chi_orb_M3_den)
    
    print(f"chi_orb(M_3) = -({B_2g}) / (4 * {g}) = {chi_orb_M3}\n")

    # Step 2: Calculate chi_orb(M_3^h)
    print("--- Step 2: Compute chi_orb(M_3^h) ---")
    print("The hyperelliptic locus M_3^h has an orbifold Euler characteristic related to that")
    print("of the moduli space of 2g+2 = 8 distinct ordered points on the projective line, M_{0,8}.")
    
    n = 2 * g + 2
    chi_M0n_val = (-1)**(n - 3) * math.factorial(n - 3)
    
    print(f"The formula is chi(M_{{0,n}}) = (-1)^(n-3) * (n-3)!.")
    print(f"For n = {n}, chi(M_{{0,8}}) = (-1)^({n-3}) * ({n-3})! = {chi_M0n_val}.")

    print("\nThe formula relating this to the hyperelliptic locus is chi_orb(M_g^h) = (1/2) * chi(M_{0,2g+2}).")
    chi_orb_M3h = Fraction(chi_M0n_val, 2)
    print(f"chi_orb(M_3^h) = (1/2) * ({chi_M0n_val}) = {chi_orb_M3h}\n")
    
    # Step 3: Final Calculation
    print("--- Step 3: Final Calculation ---")
    print("chi_orb([U/G]) = chi_orb(M_3^nh) = chi_orb(M_3) - chi_orb(M_3^h)")
    
    result = chi_orb_M3 - chi_orb_M3h
    
    print(f"\nThe final calculation is:")
    print(f"{chi_orb_M3} - ({chi_orb_M3h}) = {result}")

    return result

if __name__ == '__main__':
    final_answer = solve_euler_characteristic()
    print(f"\n<<< {final_answer.numerator}/{final_answer.denominator} >>>")