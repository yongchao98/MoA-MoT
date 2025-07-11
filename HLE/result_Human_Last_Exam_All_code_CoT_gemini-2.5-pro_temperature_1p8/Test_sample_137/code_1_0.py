import sympy
import math
from sympy.combinatorics.partitions import Partition
from sympy.combinatorics.named_groups import SymmetricGroup

def solve_hurwitz_number():
    """
    Calculates the Hurwitz number for a specific ramified cover of P^1.
    """
    d = 7
    partitions_data = [
        [4, 3],
        [5, 2],
        [2, 1, 1, 1, 1, 1],
        [2, 1, 1, 1, 1, 1]
    ]
    partitions = [Partition(p) for p in partitions_data]
    C1_part, C2_part, C3_part, C4_part = partitions

    print("Problem: Calculate the Hurwitz number counting covers P^1 -> P^1 of degree 7.")
    print(f"Degree d = {d}")
    print("Ramification profiles:")
    print(f"  λ1 = {C1_part.partition}")
    print(f"  λ2 = {C2_part.partition}")
    print(f"  λ3 = {C3_part.partition}")
    print(f"  λ4 = {C4_part.partition}")
    print("\nUsing the character-theoretic formula for Hurwitz numbers:")
    print("H = (|C1|*|C2|*|C3|^2 / (d!)^2) * Σ [χ(C1)*χ(C2)*χ(C3)^2 / χ(1)^2]")
    print("-" * 50)
    print("Step 1: Calculate the sizes of the conjugacy classes.")

    C1_size = C1_part.cardinality_of_class()
    C2_size = C2_part.cardinality_of_class()
    C3_size = C3_part.cardinality_of_class()

    print(f"Size of class C1 (partition {C1_part.partition}): |C1| = 7! / (4*3) = {C1_size}")
    print(f"Size of class C2 (partition {C2_part.partition}): |C2| = 7! / (5*2) = {C2_size}")
    print(f"Size of class C3 (partition {C3_part.partition}): |C3| = 7! / (2 * 5!) = {C3_size}")

    print("\nStep 2: Calculate the pre-factor in the formula.")
    d_factorial = math.factorial(d)
    
    pre_factor_num = C1_size * C2_size * C3_size**2
    pre_factor_den = d_factorial**2
    pre_factor = sympy.Rational(pre_factor_num, pre_factor_den)

    print(f"Pre-factor = ({C1_size} * {C2_size} * {C3_size}^2) / ({d}!)**2")
    print(f"           = {pre_factor_num} / {pre_factor_den}")
    print(f"           = {pre_factor}")
    print("-" * 50)
    
    print("Step 3: Sum over the irreducible characters of S_7.")
    
    S = SymmetricGroup(d)
    char_table = S.character_table()
    
    char_partitions = list(char_table.keys())
    id_partition = Partition([1]*d)

    total_sum = 0
    sum_terms_str = []
    
    for chi_part in sorted(char_partitions, key=lambda p: str(p.partition)):
        char_values = char_table[chi_part]
        
        chi_1 = char_values[id_partition]
        chi_C1 = char_values[C1_part]
        chi_C2 = char_values[C2_part]
        chi_C3 = char_values[C3_part]
        
        term = (chi_C1 * chi_C2 * chi_C3**2) / chi_1**2
        total_sum += term
        
        chi_part_str = str(chi_part.partition).ljust(21)
        term_str = f"({chi_C1: >2}*{chi_C2: >2}*{chi_C3: >2}²)/{chi_1: >2}² = {term}"
        print(f"Character {chi_part_str}: {term_str}")
        sum_terms_str.append(str(term))

    print("\nStep 4: Calculate the total sum.")
    print("Sum of terms = " + " + ".join(sum_terms_str))
    print(f"             = {total_sum}")
    print("-" * 50)
    
    print("Step 5: Calculate the final Hurwitz number.")
    hurwitz_number = pre_factor * total_sum
    
    print(f"Hurwitz number H = (Pre-factor) * (Sum of terms)")
    print(f"                 = {pre_factor} * {total_sum}")
    print(f"                 = {hurwitz_number}")

    if not hurwitz_number.is_Integer:
        print("\nWarning: The result is not an integer. There might be a miscalculation or a conceptual error.")
    
    print("\nFinal Answer:")
    final_int_value = int(hurwitz_number)
    print(f"The Hurwitz number is {final_int_value}.")
    return final_int_value

if __name__ == '__main__':
    result = solve_hurwitz_number()
    print(f'<<<{result}>>>')
