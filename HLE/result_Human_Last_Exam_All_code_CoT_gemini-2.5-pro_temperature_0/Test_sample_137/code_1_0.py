import math
from fractions import Fraction
from sympy.combinatorics.permutations import Partition
from sympy.combinatorics.symmetric import SymmetricGroup
from sympy.combinatorics.group_numbers import conjugacy_class_size

def solve_hurwitz_number():
    """
    Calculates the Hurwitz number for a degree 7 cover of P^1 with four specified branch points.
    """
    d = 7
    k = 4
    
    # Ramification profiles as partitions of d=7
    p1 = Partition([4, 3])
    p2 = Partition([5, 2])
    p3 = Partition([2, 1, 1, 1, 1, 1])
    partitions = [p1, p2, p3, p3]

    # The automorphism group size for the list of partitions.
    # Two partitions are the same, so it's 2! = 2.
    aut_factor = 2
    
    # Calculate sizes of the corresponding conjugacy classes in S_7
    c1_size = conjugacy_class_size(p1)
    c2_size = conjugacy_class_size(p2)
    c3_size = conjugacy_class_size(p3)
    
    d_factorial = math.factorial(d)
    
    # Get the character table of the symmetric group S_7
    S7 = SymmetricGroup(d)
    char_table = S7.character_table()
    
    # Find the column indices for the required conjugacy classes
    class_indices = {
        p: char_table.conjugacy_classes.index(p) for p in set(partitions)
    }
    id_partition = Partition([1]*d)
    id_class_idx = char_table.conjugacy_classes.index(id_partition)

    # Sum over all irreducible characters of S_7
    total_sum = 0
    for i in range(len(char_table.conjugacy_classes)):
        char_values = char_table[i]
        
        chi_id = char_values[id_class_idx]
        chi_p1 = char_values[class_indices[p1]]
        chi_p2 = char_values[class_indices[p2]]
        chi_p3 = char_values[class_indices[p3]]
        
        # The formula term for this character
        numerator = chi_p1 * chi_p2 * chi_p3**2
        denominator = chi_id**(k - 2)
        
        total_sum += numerator / denominator

    # Calculate the pre-factor for the formula
    pre_factor_numerator = c1_size * c2_size * c3_size**2
    pre_factor = pre_factor_numerator / d_factorial
    
    # Calculate the final Hurwitz number
    hurwitz_number = (1 / aut_factor) * pre_factor * total_sum
    
    # The result must be an integer
    final_result = round(hurwitz_number)

    # Print the step-by-step calculation as requested
    print("The Hurwitz number H is calculated using the formula:")
    print("H = (1/|Aut|) * (|C(4,3)| * |C(5,2)| * |C(2,1^5)|^2 / 7!) * Sum")
    print("\nWhere Sum = sum over all irreducible characters chi of S7 of:")
    print("  (chi(C(4,3)) * chi(C(5,2)) * chi(C(2,1^5))^2) / chi(id)^2\n")
    
    # Use Fraction for a more precise representation of the sum
    total_sum_fraction = Fraction(total_sum).limit_denominator()

    print("Plugging in the numbers:")
    print(f"H = (1/{aut_factor}) * ({c1_size} * {c2_size} * {c3_size}**2 / {d_factorial}) * {total_sum_fraction}")
    print(f"H = {1/aut_factor} * ({pre_factor_numerator} / {d_factorial}) * {total_sum_fraction}")
    print(f"H = {1/aut_factor} * {pre_factor} * {total_sum_fraction}")
    print(f"H = {final_result}")

solve_hurwitz_number()