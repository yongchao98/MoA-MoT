import sympy
from sympy import Rational
from math import factorial

def solve_hurwitz_number():
    """
    Calculates the Hurwitz number for a ramified cover of P^1 by P^1
    with degree 7 and four specified ramification profiles.
    """
    # 1. Define the parameters of the problem
    # Degree of the cover
    d = 7
    # Ramification profiles as partitions
    partitions_defs = [[4, 3], [5, 2], [2, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1, 1]]
    partitions = [sympy.combinatorics.partitions.Partition(p) for p in partitions_defs]
    # Number of branch points
    k = len(partitions)

    # 2. Set up the required group theory objects from sympy
    # Symmetric group S_d
    S = sympy.combinatorics.named_groups.SymmetricGroup(d)
    # Get the character table of S_d
    char_table = S.character_table()
    
    # Get the list of conjugacy classes in the same order as the table columns
    table_conjugacy_classes = char_table.conjugacy_classes
    table_partitions = [c.partition for c in table_conjugacy_classes]

    # Find the indices of the columns corresponding to our partitions
    try:
        partition_indices = [table_partitions.index(p) for p in partitions]
    except ValueError as e:
        print(f"Error: A partition was not found in the conjugacy classes of S_{d}.")
        print(e)
        return

    # Get the conjugacy class objects for our partitions
    our_conjugacy_classes = [table_conjugacy_classes[i] for i in partition_indices]
    
    # 3. Calculate the components for the Hurwitz formula
    # Calculate the sizes of these conjugacy classes
    class_sizes = [c.size for c in our_conjugacy_classes]

    # Calculate the character sum part of the formula
    char_sum = Rational(0)
    num_characters = len(char_table)

    for char_index in range(num_characters):
        # Get the degree of the character, f^chi = chi(id)
        # The identity element corresponds to the first conjugacy class (column 0).
        f_chi = char_table[char_index][0]

        # Get the character values for each of our partitions
        char_values = [char_table[char_index][j] for j in partition_indices]
        
        # Calculate the product of the character values
        prod_char_values = 1
        for val in char_values:
            prod_char_values *= val
            
        # The term in the sum is (prod chi(mu_i)) / (f^chi)^(k-2)
        term = Rational(prod_char_values, f_chi ** (k - 2))
        
        char_sum += term

    # Calculate the pre-factor (prod |C_mui|) / d!
    d_factorial = factorial(d)
    prod_class_sizes = 1
    for size in class_sizes:
        prod_class_sizes *= size
    prefactor = Rational(prod_class_sizes, d_factorial)

    # Calculate the final Hurwitz number
    hurwitz_number = prefactor * char_sum
    final_answer = int(hurwitz_number)

    # 4. Print the calculation steps and the final answer
    print(f"The degree of the cover is d = {d}.")
    print(f"The number of branch points is k = {k}.")
    print(f"The ramification profiles are {partitions_defs[0]}, {partitions_defs[1]}, {partitions_defs[2]}, and {partitions_defs[3]}.")
    print(f"The sizes of the corresponding conjugacy classes in S_{d} are {class_sizes[0]}, {class_sizes[1]}, {class_sizes[2]}, and {class_sizes[3]}.")
    print(f"d! = {d_factorial}.")
    print(f"The sum over the characters evaluates to {char_sum}.")
    
    print("\nThe final equation is:")
    print(f"H = ( {class_sizes[0]} * {class_sizes[1]} * {class_sizes[2]} * {class_sizes[3]} / {d_factorial} ) * ({char_sum}) = {final_answer}")

solve_hurwitz_number()