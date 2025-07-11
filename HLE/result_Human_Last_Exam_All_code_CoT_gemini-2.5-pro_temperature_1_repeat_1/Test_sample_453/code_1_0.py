from fractions import Fraction

def calculate_relativistic_correction():
    """
    Calculates the second-order relativistic energy correction for a hydrogen atom
    in the state n=3, l=2.
    """
    # Given quantum numbers
    n = 3
    l = 2

    # L is a shorthand for l + 1/2
    L = Fraction(l) + Fraction(1, 2)

    print("This script calculates the second-order energy shift due to the relativistic kinetic energy correction.")
    print("The calculation is based on the series expansion of the exact Klein-Gordon energy eigenvalue.")
    print("The formula for the shift (the O(α^6) term) is:")
    print("ΔE = - (m_e * c^2 * α^6 / (2 * n^3)) * [1/(4*L^3) + 3/(4*n*L^2) - 3/(2*n^2*L) + 5/(8*n^3)]")
    print("where m_e is the electron mass, c is the speed of light, and α is the fine-structure constant.")
    print("-" * 50)

    # Step-by-step calculation
    print("Step-by-step calculation for n=3, l=2:")
    print(f"1. Define quantum numbers:")
    print(f"   n = {n}")
    print(f"   l = {l}")
    print(f"   L = l + 1/2 = {l} + 1/2 = {L}")
    print("-" * 50)

    print("2. Calculate the terms inside the square bracket:")
    # Calculate each term in the bracket using exact fraction arithmetic
    term1 = Fraction(1) / (4 * L**3)
    term2 = Fraction(3) / (4 * n * L**2)
    term3 = -Fraction(3) / (2 * n**2 * L)
    term4 = Fraction(5) / (8 * n**3)

    print(f"   Term 1: 1/(4*L^3) = 1/(4*({L})^3) = {term1}")
    print(f"   Term 2: 3/(4*n*L^2) = 3/(4*{n}*({L})^2) = {term2}")
    print(f"   Term 3: -3/(2*n^2*L) = -3/(2*{n**2}*({L})) = {term3}")
    print(f"   Term 4: 5/(8*n^3) = 5/(8*{n**3}) = {term4}")
    print("-" * 50)

    print("3. Sum the terms in the bracket:")
    bracket_sum = term1 + term2 + term3 + term4
    print(f"   Sum = {term1} + {term2} + {term3} + {term4} = {bracket_sum}")
    print("-" * 50)

    print("4. Calculate the prefactor outside the bracket:")
    prefactor = -Fraction(1) / (2 * n**3)
    print(f"   Prefactor = -1/(2*n^3) = -1/(2*{n**3}) = {prefactor}")
    print("-" * 50)

    print("5. Calculate the final coefficient for the energy shift:")
    final_coeff = prefactor * bracket_sum
    print(f"   Total Coefficient = (Prefactor) * (Bracket Sum)")
    print(f"                     = ({prefactor}) * ({bracket_sum}) = {final_coeff}")
    print("-" * 50)

    print("\nThe final expression for the energy shift is:")
    print(f"ΔE = ({final_coeff.numerator}/{final_coeff.denominator}) * m_e * c^2 * α^6")


if __name__ == "__main__":
    calculate_relativistic_correction()
    # The final expression is the answer.
    # For n=3, l=2, L=5/2.
    # term1 = 2/125, term2 = 1/25, term3 = -1/15, term4 = 5/216.
    # bracket_sum = -4/375 + 5/216 = 337/27000.
    # prefactor = -1/54.
    # final_coeff = (-1/54) * (337/27000) = -337/1458000.
    final_expression = "(-337/1458000) * m_e * c^2 * α^6"
    # print(f"\n<<<{final_expression}>>>") # Suppressed for final output format