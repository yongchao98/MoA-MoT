from fractions import Fraction

def calculate_relativistic_energy_shift():
    """
    Calculates the second-order relativistic energy shift for a hydrogen atom
    in the state n=3, l=2.

    The calculation is based on the difference between the alpha^6 term from the
    Klein-Gordon equation's exact solution and the expectation value of the p^6
    term of the kinetic energy expansion.
    
    ΔE_n^(2) = ΔE(α^6) - <n|H''|n>
    """
    n = 3
    l = 2

    # Use Fractions for exact arithmetic
    n_f = Fraction(n)
    l_f = Fraction(l)
    l_prime = l_f + Fraction(1, 2) # l + 1/2

    # Part 1: Calculate the coefficient of ΔE(α^6)
    # ΔE(α^6) = -mc^2 α^6 * [ 1/(8n^3(l+1/2)^3) + 3/(8n^4(l+1/2)^2) 
    #                       - 3/(4n^5(l+1/2)) + 5/(16n^6) ]
    term1_E6 = Fraction(1) / (8 * n_f**3 * l_prime**3)
    term2_E6 = Fraction(3) / (8 * n_f**4 * l_prime**2)
    term3_E6 = -Fraction(3) / (4 * n_f**5 * l_prime)
    term4_E6 = Fraction(5) / (16 * n_f**6)
    
    coeff_E6 = -(term1_E6 + term2_E6 + term3_E6 + term4_E6)

    # Part 2: Calculate the coefficient of <n|H''|n>
    # <n|H''|n> = (mc^2 α^6 / 2) * [ 5/(8n^6) - 3/(2n^5(l+1/2)) 
    #                              - 1/(n^3*l*(l+1/2)*(l+1)) ]
    term1_Hpp = Fraction(5) / (8 * n_f**6)
    term2_Hpp = -Fraction(3) / (2 * n_f**5 * l_prime)
    term3_Hpp = -Fraction(1) / (n_f**3 * l_f * l_prime * (l_f + 1))
    
    coeff_Hpp = Fraction(1, 2) * (term1_Hpp + term2_Hpp + term3_Hpp)

    # Part 3: Calculate the final coefficient for ΔE_n^(2)
    final_coeff = coeff_E6 - coeff_Hpp
    
    num = final_coeff.numerator
    den = final_coeff.denominator
    c_exp = 2
    a_exp = 6

    print("The second-order energy shift due to the relativistic correction is given by the equation:")
    print(f"ΔE = ({num} / {den}) * m * c^{c_exp} * α^{a_exp}")
    print("\nWhere:")
    print("m = mass of the electron")
    print("c = speed of light")
    print("α = fine-structure constant")

calculate_relativistic_energy_shift()
<<<ΔE = (1319 / 729000) * m * c^2 * α^6>>>