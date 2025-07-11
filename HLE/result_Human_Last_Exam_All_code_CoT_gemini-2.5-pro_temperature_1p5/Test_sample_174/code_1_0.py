from sympy import bernoulli
from fractions import Fraction

def chi_orb_Mg(g):
    """
    Computes the orbifold Euler characteristic of the moduli stack of
    genus g curves, M_g, using the Harer-Zagier formula.
    Formula: chi(M_g) = B_{2g} / (2g * (2g - 2)) for g > 1.
    """
    if g <= 1:
        raise ValueError("Formula is for g > 1")
    # The bernoulli function in sympy can return fractions, so we create a Fraction object.
    b_2g = Fraction(bernoulli(2*g))
    denominator = 2*g * (2*g - 2)
    return b_2g / denominator

def chi_orb_Hg(g):
    """
    Computes the orbifold Euler characteristic of the moduli stack of
    genus g hyperelliptic curves, H_g, using Getzler's formula.
    Formula: chi(H_g) = -B_{2g} / (2g) for g >= 1.
    """
    if g < 1:
        raise ValueError("Formula is for g >= 1")
    b_2g = Fraction(bernoulli(2*g))
    denominator = 2*g
    return -b_2g / denominator

# A smooth plane quartic curve has genus g=3.
g = 3

print(f"The genus of a smooth plane quartic (degree 4) is g = (4-1)(4-2)/2 = {g}.")
print("The stack of smooth plane quartics [U/G] corresponds to the stack of non-hyperelliptic genus 3 curves.")
print("Its orbifold Euler characteristic can be found by subtracting the characteristic of the hyperelliptic locus H_3")
print("from the characteristic of the entire moduli stack of genus 3 curves M_3.")
print("-" * 20)

# Step 1: Calculate chi_orb(M_3)
chi_M3 = chi_orb_Mg(g)
b6_val = bernoulli(2*g)
b6_frac = Fraction(b6_val)
print(f"Step 1: Compute chi(M_3) using the Harer-Zagier formula chi(M_g) = B_{2*g} / (2g * (2g-2)).")
print(f"For g=3, B_6 = {b6_frac}. So, chi(M_3) = ({b6_frac}) / (6 * 4) = {chi_M3}.")
print("-" * 20)

# Step 2: Calculate chi_orb(H_3)
chi_H3 = chi_orb_Hg(g)
print(f"Step 2: Compute chi(H_3) using Getzler's formula chi(H_g) = -B_{2*g} / (2g).")
print(f"For g=3, chi(H_3) = -({b6_frac}) / 6 = {chi_H3}.")
print("-" * 20)

# Step 3: Calculate the final result using the additivity principle.
result = chi_M3 - chi_H3
print("Step 3: Compute the final result.")
print("chi([U/G]) = chi(M_3) - chi(H_3)")
print(f"           = {chi_M3} - ({chi_H3})")
print(f"           = {chi_M3} + {-chi_H3}")
# To show the common denominator calculation
common_den = result.denominator
term1_adj = chi_M3.numerator * (common_den // chi_M3.denominator)
term2_adj = (-chi_H3).numerator * (common_den // (-chi_H3).denominator)
print(f"           = {term1_adj}/{common_den} + {term2_adj}/{common_den}")
print(f"           = {result}")
