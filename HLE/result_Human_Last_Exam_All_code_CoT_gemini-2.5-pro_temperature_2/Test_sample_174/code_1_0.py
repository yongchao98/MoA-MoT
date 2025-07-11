from fractions import Fraction

def compute_euler_characteristic():
    """
    Computes the orbifold Euler characteristic of the moduli stack of smooth plane quartics.
    """
    g = 3
    print(f"The genus of a smooth plane quartic is g = {g}.")
    print("The stack of smooth plane quartics is the non-hyperelliptic locus in the moduli stack of genus 3 curves, M_3.")
    print("We use the formula: chi_orb([U/G]) = chi_orb(M_3) - chi_orb(H_3), where H_3 is the hyperelliptic locus.\n")

    # The 6th Bernoulli number, B_6
    B6 = Fraction(1, 42)
    print(f"The 6th Bernoulli number B_6 is {B6}.")

    # Calculate chi_orb(M_3) = zeta(1-2g) = zeta(-5) = -B_6 / 6
    zeta_minus_5 = -B6 / 6
    chi_orb_M3 = zeta_minus_5
    print(f"The orbifold Euler characteristic of M_3 is zeta(-5) = -B_6/6 = {chi_orb_M3}.")

    # Calculate chi_orb(H_3) = zeta(-5) / (1 - 2^(2g-1))
    denominator = 1 - 2**(2 * g - 1)
    chi_orb_H3 = zeta_minus_5 / denominator
    print(f"The orbifold Euler characteristic of the hyperelliptic locus H_3 is zeta(-5)/(1-2^(2g-1)) = {zeta_minus_5} / {denominator} = {chi_orb_H3}.\n")
    
    # Calculate the final result
    result = chi_orb_M3 - chi_orb_H3

    print(f"The final result is the difference between these two values.")
    print(f"chi_orb([U/G]) = chi_orb(M_3) - chi_orb(H_3)")
    print(f"                = {chi_orb_M3} - ({chi_orb_H3})")
    print(f"                = {result}")

compute_euler_characteristic()