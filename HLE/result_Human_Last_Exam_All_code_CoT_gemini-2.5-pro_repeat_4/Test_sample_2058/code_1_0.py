from fractions import Fraction

def solve_mass_formula():
    """
    Calculates the total mass M(A_5, rho, 2) as defined in the problem.
    """
    
    # The group is G = A_5, its order is |A_5| = 60.
    # The prime is p = 2.
    # The representation is the 5-dim permutation representation rho, which decomposes
    # as rho = 1 + rho_4, where 1 is the trivial representation and rho_4 is a
    # 4-dimensional irreducible representation.
    # The conductor c(rho . phi) is equal to c(rho_4 . phi), since c(1 . phi) = 0.
    
    # We will compute the total sum S = sum_{phi: Gamma_2 -> A_5} 2**(-c(rho_4 . phi)).
    # The sum is partitioned by the image of the homomorphism phi.
    
    # Initialize the total sum S.
    total_sum_S = Fraction(0)
    
    # 1. Contribution from unramified homomorphisms.
    # An unramified homomorphism phi is determined by the image of Frob_2.
    # There are |A_5| = 60 such homomorphisms. For all of them, the conductor is 0.
    # Their contribution to the sum S is 60 * 2**0 = 60.
    unramified_contribution = Fraction(60)
    total_sum_S += unramified_contribution
    
    # 2. Contribution from ramified homomorphisms.
    # The image H of a ramified homomorphism from Gamma_2 must have a normal
    # 2-subgroup K such that H/K is cyclic. For subgroups of A_5, the only
    # possible images for ramified maps are C_2, V_4, and A_4.
    
    # 2a. Image H = C_2 (cyclic group of order 2).
    # There are 15 subgroups of A_5 isomorphic to C_2.
    # There are 6 ramified quadratic extensions of Q_2, which correspond to
    # 6 ramified epimorphisms epsilon: Gamma_2 -> C_2.
    # The conductor of the composed representation rho_4 . phi is 2*f(epsilon),
    # where f(epsilon) is the conductor of the quadratic character epsilon.
    # - 2 of these extensions have f=2. Contribution: 15 * 2 * 2**(-2*2) = 30/16 = 15/8.
    total_sum_S += Fraction(15, 8)
    # - 4 of these extensions have f=3. Contribution: 15 * 4 * 2**(-2*3) = 60/64 = 15/16.
    total_sum_S += Fraction(15, 16)

    # 2b. Image H = V_4 (Klein four-group).
    # There are 30 injective homomorphisms from V_4 to A_5.
    # There are 3 distinct V_4 extensions of Q_2. For each extension, we have
    # a set of corresponding maps to A_5.
    # The conductor c(rho_4 . phi) is the sum of the conductors of the three
    # quadratic characters associated with the subfields of the V_4 extension.
    # - Extension 1: Q_2(sqrt(3), sqrt(5)). Conductors of subfields (2, 0, 2). Total c = 4.
    #   Contribution: 30 * 2**(-4) = 30/16 = 15/8.
    total_sum_S += Fraction(15, 8)
    # - Extension 2: Q_2(sqrt(2), sqrt(5)). Conductors (3, 0, 3). Total c = 6.
    #   Contribution: 30 * 2**(-6) = 30/64 = 15/32.
    total_sum_S += Fraction(15, 32)
    # - Extension 3: Q_2(sqrt(2), sqrt(-1)). Conductors (3, 2, 3). Total c = 8.
    #   Contribution: 30 * 2**(-8) = 30/256 = 15/128.
    total_sum_S += Fraction(15, 128)

    # 2c. Image H = A_4.
    # There are 5 subgroups of A_5 isomorphic to A_4. The number of epimorphisms
    # from Gamma_2 to A_4 is 24. This gives 5 * 24 = 120 maps with image A_4.
    # There is one A_4 extension of Q_2. The conductor of the relevant 3D part of
    # the representation is 4. So c(rho_4 . phi) = 4.
    # Contribution: 120 * 2**(-4) = 120/16 = 15/2.
    total_sum_S += Fraction(15, 2)

    # The order of the group A_5 is 60.
    G_order = 60
    
    # The total mass is M = S / |A_5|.
    M = total_sum_S / G_order
    
    # Print the final result.
    print("The calculation for the total mass M(A_5, rho, 2) results in a fraction.")
    print(f"The numerator of the fraction is: {M.numerator}")
    print(f"The denominator of the fraction is: {M.denominator}")
    print(f"The final equation is: M(A_5, rho, 2) = {M.numerator} / {M.denominator}")

solve_mass_formula()
<<<621/512>>>