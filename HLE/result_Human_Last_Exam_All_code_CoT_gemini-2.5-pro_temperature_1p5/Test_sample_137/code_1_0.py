import math
from fractions import Fraction
from collections import Counter

def conjugacy_class_size(partition, d):
    """Calculates the size of a conjugacy class in S_d."""
    counts = Counter(partition)
    denominator = 1
    for k, m_k in counts.items():
        denominator *= (k**m_k) * math.factorial(m_k)
    return math.factorial(d) // denominator

def solve_hurwitz_number():
    """
    Calculates the Hurwitz number H((4,3), (5,2), (2,1^5), (2,1^5)) for degree 7 covers.
    """
    d = 7
    p1 = (4, 3)
    p2 = (5, 2)
    p3 = (2, 1, 1, 1, 1, 1)
    
    # In ATLAS of Finite Groups, these classes for S_7 are:
    # (4,3) -> 4B
    # (5,2) -> 5B
    # (2,1^5) -> 2A
    
    # Character table data for S_7 from ATLAS.
    # Each tuple: (dimension, chi(2A), chi(4B), chi(5B))
    # which corresponds to (chi(id), chi(2,1^5), chi(4,3), chi(5,2))
    character_data = [
        (1, 1, 1, 1),
        (6, 4, 0, -1),
        (14, 6, 0, -1),
        (14, -2, 2, 1),
        (15, 7, -1, 0),
        (15, -1, -1, 0),
        (20, -4, 0, 0),
        (21, 5, 1, 1),
        (21, -3, 1, 1),
        (35, 7, -1, 0),
        (35, -5, -1, 0),
        (1, -1, -1, -1), # Sign character
        (6, -4, 0, 1),
        (14, 2, -2, -1),
        (20, 4, 0, 0)
    ]

    # Calculate sizes of conjugacy classes
    size_p1 = conjugacy_class_size(p1, d)
    size_p2 = conjugacy_class_size(p2, d)
    size_p3 = conjugacy_class_size(p3, d)

    d_factorial = math.factorial(d)
    
    # Calculate the pre-factor in the Hurwitz number formula
    prefactor = (size_p1 * size_p2 * size_p3 * size_p3) // d_factorial
    
    # Calculate the sum over characters
    character_sum = Fraction(0)
    sum_terms = []
    
    for i, (dim, chi_p3, chi_p1, chi_p2) in enumerate(character_data):
        if dim == 0: continue
        
        numerator = chi_p1 * chi_p2 * chi_p3 * chi_p3
        denominator = dim * dim
        
        term = Fraction(numerator, denominator)
        character_sum += term
        sum_terms.append(f"({term})")

    # The Hurwitz number is the product
    hurwitz_number = prefactor * character_sum
    
    # Print the explanation and the final result
    print("This program calculates the Hurwitz number counting ramified covers P^1 -> P^1 of degree 7")
    print("with four branch points having ramification profiles (4,3), (5,2), (2,1^5), (2,1^5).\n")
    print("The formula used is H = ( |C1|*|C2|*|C3|*|C4| / d! ) * Sum_chi [ chi(C1)*chi(C2)*chi(C3)*chi(C4) / chi(id)^2 ]\n")
    print(f"Degree d = {d}")
    print(f"Partitions: lambda_1={p1}, lambda_2={p2}, lambda_3={p3}, lambda_4={p3}\n")
    
    print("Sizes of the conjugacy classes in S_7:")
    print(f"|C(4,3)| = 7! / (4^1 * 1! * 3^1 * 1!) = {size_p1}")
    print(f"|C(5,2)| = 7! / (5^1 * 1! * 2^1 * 1!) = {size_p2}")
    print(f"|C(2,1^5)| = 7! / (2^1 * 1! * 1^5 * 5!) = {size_p3}\n")
    
    print("Calculation:")
    print(f"Prefactor = (|C(4,3)| * |C(5,2)| * |C(2,1^5)|^2) / 7!")
    print(f"Prefactor = ({size_p1} * {size_p2} * {size_p3}^2) / {d_factorial} = {prefactor}\n")
    
    sum_str = " + ".join(sum_terms)
    # The string for the sum can be very long. Let's show only non-zero terms for clarity.
    non_zero_sum_terms = [s for s in sum_terms if s != '(0)']
    print(f"Sum over characters S = {' + '.join(non_zero_sum_terms)}")
    print(f"S = {character_sum}\n")
    
    print(f"Hurwitz Number H = Prefactor * S")
    print(f"H = {prefactor} * {character_sum} = {hurwitz_number}")
    print("\nThe final result is:")
    print(int(hurwitz_number))

solve_hurwitz_number()
<<<39984>>>