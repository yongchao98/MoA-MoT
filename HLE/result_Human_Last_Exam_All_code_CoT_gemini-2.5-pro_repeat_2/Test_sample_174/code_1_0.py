from fractions import Fraction

def solve_euler_characteristic():
    """
    Calculates the orbifold Euler characteristic of the moduli stack of smooth plane quartics.
    """
    # The genus g of a smooth plane quartic curve (degree d=4) is g = (d-1)(d-2)/2.
    g = 3
    
    # The relevant Bernoulli number is B_{2g} = B_6.
    # The value of the 6th Bernoulli number is B_6 = 1/42.
    B6 = Fraction(1, 42)
    
    # The orbifold Euler characteristic of the moduli stack of genus g curves, M_g,
    # is given by the Harer-Zagier formula: chi(M_g) = -B_{2g} / (2g).
    chi_M3 = -B6 / (2 * g)
    
    # The orbifold Euler characteristic of the hyperelliptic locus, H_g,
    # is given by the formula: chi(H_g) = (g-1) * B_{2g} / (2g).
    chi_H3 = (g - 1) * B6 / (2 * g)
    
    # The moduli stack of smooth plane quartics [U/G] is the moduli stack of
    # non-hyperelliptic curves of genus 3, M_3^nh.
    # Its Euler characteristic is chi(M_3) - chi(H_3).
    result = chi_M3 - chi_H3

    print("The problem asks for the orbifold Euler characteristic of the quotient stack [U/G] of smooth plane quartics.")
    print("This space is the moduli stack of non-hyperelliptic curves of genus 3, denoted M_3^nh.")
    print("Its Euler characteristic is calculated as the difference between the Euler characteristic of the full moduli stack M_3 and the hyperelliptic locus H_3.")
    print("\nCalculation steps:")
    print(f"The genus g of a smooth plane quartic is g = (4-1)(4-2)/2 = {g}.")
    print(f"The required Bernoulli number is B_{2*g} = B_6, which is {B6}.")
    print(f"The orbifold Euler characteristic of M_3 is: χ(M_3) = -B_6 / 6 = -({B6}) / 6 = {chi_M3}.")
    print(f"The orbifold Euler characteristic of H_3 is: χ(H_3) = (3-1)*B_6 / 6 = 2 * ({B6}) / 6 = {chi_H3}.")
    print(f"\nThe final result is: χ([U/G]) = χ(M_3) - χ(H_3) = {chi_M3} - ({chi_H3}) = {result}.")

solve_euler_characteristic()