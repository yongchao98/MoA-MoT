import sympy
from sympy.combinatorics.partitions import Partition
from fractions import Fraction

def solve_hurwitz_number():
    """
    Calculates the Hurwitz number for a specific ramified cover of P^1.

    The problem is to find the number of degree 7 covers of P^1 with four branch points
    having ramification profiles (4,3), (5,2), (2,1,1,1,1,1), and (2,1,1,1,1,1).
    This number corresponds to the count of monodromy tuples in S_7.
    """
    # 1. Define parameters of the problem
    d = 7
    partitions = [
        Partition([4, 3]),
        Partition([5, 2]),
        Partition([2, 1, 1, 1, 1, 1]),
        Partition([2, 1, 1, 1, 1, 1])
    ]
    k = len(partitions)
    d_factorial = sympy.factorial(d)

    # 2. Get the character table and conjugacy classes for S_7
    S_d = sympy.symmetrics.SymmetricGroup(d)
    char_table = CharacterTable(S_d)

    # 3. Find class sizes and indices for the given partitions
    class_indices = []
    class_sizes = []
    class_map = {c.partition: c for c in char_table.conjugacy_classes}
    for p in partitions:
        c_class = class_map[p]
        class_indices.append(char_table.conjugacy_classes.index(c_class))
        class_sizes.append(len(c_class))
    
    # 4. Compute the sum over all irreducible characters
    sum_of_terms = Fraction(0)
    for char in char_table:
        chi_id = char[0]  # Dimension of the character
        
        # Numerator is the product of character values on the specified classes
        numerator = 1
        for idx in class_indices:
            numerator *= char[idx]
        
        # Denominator depends on the character dimension and number of branch points
        denominator = chi_id**(k - 2)
        
        if denominator != 0:
            term = Fraction(numerator, denominator)
            sum_of_terms += term

    # 5. Calculate the total number of tuples (N)
    prod_class_sizes = 1
    for size in class_sizes:
        prod_class_sizes *= size

    # The formula for N is (prod(|C_i|) / d!) * sum_of_terms
    N = Fraction(prod_class_sizes, d_factorial) * sum_of_terms

    # 6. Print the result in a readable equation format
    # The final Hurwitz number is an integer, so we expect the fraction to resolve
    print(f"The number of covers is calculated using the formula:")
    print(f"N = (|C_(4,3)| * |C_(5,2)| * |C_(2,1^5)|^2) / 7! * Sum_of_Terms")
    print(f"N = ({class_sizes[0]} * {class_sizes[1]} * {class_sizes[2]}^2) / {d_factorial} * ({sum_of_terms})")
    print(f"N = {int(N)}")
    
solve_hurwitz_number()