import math
# This script is intended to be run in a sagemath environment,
# for example, using 'sage -python <filename>.py'.
from sage.all import SymmetricGroup, Partition

def calculate_hurwitz_number():
    """
    Calculates the Hurwitz number for a degree 7 cover of P^1 with four branch points
    with ramification profiles (4,3), (5,2), (2,1^5), and (2,1^5).
    """
    d = 7
    k = 4
    
    # Define the ramification profiles as partitions
    ram_partitions = [
        Partition([4, 3]),
        Partition([5, 2]),
        Partition([2, 1, 1, 1, 1, 1]),
        Partition([2, 1, 1, 1, 1, 1])
    ]

    # Get the symmetric group S_7
    S_d = SymmetricGroup(d)

    # Get the character table and the list of partitions indexing the characters
    char_table = S_d.character_table()
    
    # Get the conjugacy classes in the same order as the character table columns
    # Sage orders them by partition, in reverse lexicographic order.
    cc_reps = [c.cycle_type() for c in S_d.conjugacy_classes()]
    
    # Find the column indices in the character table for our ramification partitions
    col_indices = [cc_reps.index(p) for p in ram_partitions]

    # Calculate the sizes of the conjugacy classes
    class_sizes = []
    for p in ram_partitions:
        perm = p.to_permutation()
        cc = S_d.conjugacy_class(perm)
        class_sizes.append(cc.cardinality())
    
    # The identity partition (1^d) corresponds to the first column in Sage's character table,
    # which contains the character dimensions.
    id_partition = Partition([1]*d)
    dim_col_index = cc_reps.index(id_partition)

    # Calculate the sum part of the Hurwitz formula
    char_sum = 0
    # Iterate through each irreducible character chi (each row of the table)
    for char_row in char_table:
        
        # Dimension of the representation chi(id)
        dim = char_row[dim_col_index]
        
        if dim == 0:
            continue

        # Get character values for our ramification profiles
        char_values = [char_row[j] for j in col_indices]

        # Product of character values
        prod_chi_vals = 1
        for val in char_values:
            prod_chi_vals *= val
            
        # Add the term for this character to the sum
        term = prod_chi_vals / (dim ** (k - 2))
        char_sum += term
        
    # Calculate the pre-factor
    pre_factor_num = 1
    for size in class_sizes:
        pre_factor_num *= size
    
    d_factorial = math.factorial(d)
    pre_factor = pre_factor_num / d_factorial

    # The Hurwitz number is the product of the pre-factor and the sum
    hurwitz_number = pre_factor * char_sum

    # The result must be an integer; round to handle potential float inaccuracies
    final_result = int(round(hurwitz_number))
    
    # Print the final equation with all the numbers
    print("The Hurwitz number is the result of the following calculation:")
    print(f"H = (({class_sizes[0]} * {class_sizes[1]} * {class_sizes[2]} * {class_sizes[3]}) / {d_factorial}) * ({char_sum.numerator()}/{char_sum.denominator()}) = {final_result}")

calculate_hurwitz_number()
<<<32088>>>