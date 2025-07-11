import math
from fractions import Fraction

def calculate_total_mass():
    """
    Calculates the total mass M(G, rho, p) for G=A_5, rho the 5-dimensional 
    permutation representation, and p=2, based on Kedlaya's mass formula.

    The formula is:
    M(G, rho, p) = sum_{g in G, p does not divide ord(g)} (m_p(<g>)/|Z_G(g)|) * p^(-a(rho|<g>))

    For this problem (G=A_5, p=2, rho=permutation rep), this simplifies to:
    M(A_5, rho, 2) = sum_{g in A_5, ord(g) is odd} (1 / |Z_{A_5}(g)|) * 2^-(5 - chi_rho(g))

    The sum is computed by iterating through the conjugacy classes of A_5 with odd order.
    """
    
    # Data for the conjugacy classes of A_5 with odd order elements
    # Structure: (class_name, num_elements, centralizer_size, char_rho)
    # char_rho is the number of fixed points for the permutation representation.
    classes_data = [
        ('identity', 1, 60, 5),
        ('3-cycles', 20, 3, 2),
        ('5-cycles', 24, 5, 0)
    ]

    total_mass = Fraction(0)
    
    print("Calculating the total mass M(A_5, rho, 2) using Kedlaya's formula.")
    print("The total mass is the sum of contributions from all elements of odd order.")
    print("We group these elements by conjugacy class.\n")

    dim_rho = 5
    p = 2
    equation_parts = []

    for name, num_elements, centralizer_size, char_rho in classes_data:
        conductor = dim_rho - char_rho
        
        # Use Fraction for exact arithmetic
        term_per_element = Fraction(1, centralizer_size) * Fraction(1, p**conductor)
        class_contribution = num_elements * term_per_element
        
        total_mass += class_contribution
        equation_parts.append(f"({num_elements} * 1/{centralizer_size} * {p}^-({dim_rho}-{char_rho}))")
        
        print(f"Contribution from the class of {name}:")
        print(f"  Number of elements: {num_elements}")
        print(f"  Centralizer size |Z(g)|: {centralizer_size}")
        print(f"  Character value chi(g): {char_rho}")
        print(f"  Term for one element = 1/{centralizer_size} * {p}^-({dim_rho} - {char_rho}) = {term_per_element}")
        print(f"  Total contribution for this class = {num_elements} * {term_per_element} = {class_contribution}\n")

    final_equation = f"M(A_5, rho, 2) = {' + '.join(equation_parts)}"
    # Re-calculate the numerical parts for the final equation string
    parts_values = []
    for data in classes_data:
        num, cent, char = data[1], data[2], data[3]
        contr = Fraction(num, cent * (p**(dim_rho-char)))
        parts_values.append(f"{contr.numerator}/{contr.denominator}")

    print(f"Final Calculation:")
    print(f"{' + '.join(parts_values)} = {total_mass}")
    
    # Final result in lowest terms
    final_fraction = total_mass
    print(f"\nThe total mass M(A_5, rho, 2) is {final_fraction.numerator}/{final_fraction.denominator}.")

calculate_total_mass()