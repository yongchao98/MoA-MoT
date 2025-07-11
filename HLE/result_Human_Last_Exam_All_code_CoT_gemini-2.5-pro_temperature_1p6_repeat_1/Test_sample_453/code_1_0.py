import fractions

def calculate_energy_shift():
    """
    Calculates the second-order energy shift for an electron in a hydrogen atom
    due to the relativistic kinetic energy correction.
    """
    # State quantum numbers
    n = 3
    l = 2

    print("Step 1: Calculate the first-order energy shift, ΔE_n^(1).")
    print(f"For a state (n={n}, l={l}), the first-order shift due to H' = -p^4/(8 m^3 c^2) is given by:")
    print("ΔE_n^(1) = (E_n^2 / (2*m*c^2)) * (3 - 4*n / (l + 1/2))")
    
    # Calculate the coefficient C1 = (3 - 4n/(l+1/2))/2
    # ΔE_n^(1) = C1 * E_n^2 / (m*c^2)
    frac_l_plus_half = fractions.Fraction(l * 2 + 1, 2)
    term2 = fractions.Fraction(4 * n) / frac_l_plus_half
    c1_numerator = 3 - term2
    C1 = c1_numerator / 2

    print(f"For n={n}, l={l}:")
    print(f"ΔE_3^(1) = (E_3^2 / (2*m*c^2)) * (3 - 4*{n} / ({l} + 1/2))")
    print(f"       = (E_3^2 / (2*m*c^2)) * (3 - {term2.numerator}/{term2.denominator})")
    print(f"       = (E_3^2 / (2*m*c^2)) * ({c1_numerator.numerator}/{c1_numerator.denominator})")
    print(f"       = ({C1.numerator}/{C1.denominator}) * E_3^2 / (m*c^2)")
    print("-" * 30)

    print("Step 2: Use the relation ΔE_n^(2) = - (E_n / (m*c^2)) * ΔE_n^(1) to find the second-order shift.")
    # The coefficient for the second order shift is -C1
    C2 = -C1
    
    print(f"ΔE_3^(2) = - (E_3 / (m*c^2)) * (({C1.numerator}/{C1.denominator}) * E_3^2 / (m*c^2))")
    print(f"       = ({C2.numerator}/{C2.denominator}) * E_3^3 / (m*c^2)^2")
    print("-" * 30)

    print("Step 3: Express the result in terms of fundamental constants.")
    print("Substitute the unperturbed energy E_n = -m*c^2*α^2 / (2*n^2).")
    
    # Coefficient for E_n
    E3_coeff = fractions.Fraction(-1, 2 * n**2)

    print(f"For n={n}, E_3 = ({E3_coeff.numerator}/{E3_coeff.denominator}) * m*c^2*α^2")
    
    final_coeff = C2 * (E3_coeff**3)
    
    print(f"ΔE_3^(2) = ({C2.numerator}/{C2.denominator}) * (({E3_coeff.numerator}/{E3_coeff.denominator}) * m*c^2*α^2)^3 / (m*c^2)^2")
    print(f"       = ({C2.numerator}/{C2.denominator}) * (({E3_coeff.numerator**3})/({E3_coeff.denominator**3})) * (m*c^2*α^6)")
    print(f"       = ({final_coeff.numerator}/{final_coeff.denominator}) * m*c^2*α^6")
    
    print("\nFinal Answer:")
    final_expression = f"({final_coeff.numerator}/{final_coeff.denominator}) * m * c^2 * α^6"
    print(final_expression)
    
    # Required final answer format
    print(f"\n<<<{final_expression}>>>")


calculate_energy_shift()