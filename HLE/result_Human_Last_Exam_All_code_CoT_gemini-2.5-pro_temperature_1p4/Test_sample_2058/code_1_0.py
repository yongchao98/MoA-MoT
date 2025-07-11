import math
from fractions import Fraction

def solve_mass_formula():
    """
    This function calculates the total mass M(A5, rho, 2) using Kedlaya's formula.
    M(G, rho, p) = sum_{H cyclic subgroup of G} |Hom(Z_p, H)|/|H| * p^(-c(H))
    where c(H) = dim(V) - dim(V^{H/H_p}) for the tame part of the representation.
    Here, G = A5, rho is the 5D permutation representation, and p = 2.
    """
    
    print("This script calculates the total mass M(A5, rho, 2).")
    print("The calculation is based on Kedlaya's mass formula, summing contributions from all cyclic subgroups of A5.\n")
    
    total_mass = Fraction(0)

    # Case 1: H = C_1 (trivial subgroup)
    num_H = 1
    order_H = 1
    # |Hom(Z_2, C_1)| is number of elements of 2-power order in C_1, which is 1 (the identity)
    hom_H = 1
    # c(H) = 5 - dim(V^{C_1}) = 5 - 5 = 0
    exp_c = 0
    contrib_C1 = num_H * Fraction(hom_H, order_H) * Fraction(1, 2**exp_c)
    total_mass += contrib_C1

    # Case 2: H = C_2 (subgroup of order 2)
    num_H = 15
    order_H = 2
    # |Hom(Z_2, C_2)| is number of elements of 2-power order in C_2, which is 2 (id, generator)
    hom_H = 2
    # H/H_2 = C_1, so c(H) = 5 - dim(V^{C_1}) = 0
    exp_c = 0
    contrib_C2 = num_H * Fraction(hom_H, order_H) * Fraction(1, 2**exp_c)
    total_mass += contrib_C2

    # Case 3: H = C_3 (subgroup of order 3)
    num_H = 10
    order_H = 3
    # |Hom(Z_2, C_3)| is 1 (the identity)
    hom_H = 1
    # H/H_2 = C_3. Character values for C_3 are chi(id)=5, chi(g)=2, chi(g^2)=2
    # dim(V^{C_3}) = (5+2+2)/3 = 3. So, c(H) = 5 - 3 = 2
    exp_c = 2
    contrib_C3 = num_H * Fraction(hom_H, order_H) * Fraction(1, 2**exp_c)
    total_mass += contrib_C3

    # Case 4: H = C_5 (subgroup of order 5)
    num_H = 6
    order_H = 5
    # |Hom(Z_2, C_5)| is 1 (the identity)
    hom_H = 1
    # H/H_2 = C_5. Character values are chi(id)=5, chi(g^k)=0 for k=1..4
    # dim(V^{C_5}) = (5+0+0+0+0)/5 = 1. So, c(H) = 5 - 1 = 4
    exp_c = 4
    contrib_C5 = num_H * Fraction(hom_H, order_H) * Fraction(1, 2**exp_c)
    total_mass += contrib_C5
    
    print("The total mass is the sum of contributions from each class of cyclic subgroups:")
    print(f"Contribution from {num_C1} C1 subgroup: {contrib_C1}")
    print(f"Contribution from {num_C2} C2 subgroups: {contrib_C2}")
    print(f"Contribution from {num_C3} C3 subgroups: {contrib_C3}")
    print(f"Contribution from {num_C5} C5 subgroups: {contrib_C5}")
    
    print("\nThe final equation for the total mass is:")
    # Using Fraction to display intermediate values
    val1 = contrib_C1
    val2 = contrib_C2
    val3 = contrib_C3
    val4 = contrib_C5
    
    print(f"M(A5, rho, 2) = {val1} + {val2} + {val3} + {val4}")
    print(f"M(A5, rho, 2) = {total_mass}")
    
    print(f"\nThe final result as a fraction in lowest terms is: {total_mass.numerator}/{total_mass.denominator}")

solve_mass_formula()
<<<2029/120>>>