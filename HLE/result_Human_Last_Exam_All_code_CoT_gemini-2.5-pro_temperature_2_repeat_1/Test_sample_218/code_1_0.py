import sympy
from sympy.combinatorics.partitions import Partitions

def solve_sn_irreps_dimension():
    """
    Calculates the number of irreducible representations of S_n with dimension
    below a given limit.

    This problem concerns n=25 and a limit of 500,000.
    """
    n = 25
    dimension_limit = 500000
    
    # The number of partitions of 25, p(25), is 1958.
    # We will iterate through all of them.
    partitions_generator = Partitions(n)
    
    count = 0
    print(f"Finding irreps of S_{n} with dimension < {dimension_limit}:\n")
    
    # We will print the partitions and their dimensions that meet the criteria.
    # The format will be: Partition: [p_1, p_2, ...], Dimension: D
    # This fulfills the requirement to "output each number".
    
    for p in partitions_generator:
        # The p object from sympy.combinatorics.Partitions has a dimension() method
        # which calculates the dimension of the irrep for that partition using
        # the Hook Length Formula.
        dim = p.dimension()
        
        if dim < dimension_limit:
            count += 1
            # For each partition that satisfies the condition, we print its details.
            # A partition is represented as a list of integers.
            # e.g., for n=4, a partition is [3, 1].
            partition_list = p.partition
            print(f"Partition: {partition_list}, Dimension: {dim}")

    print(f"\nTotal number of irreducible representations of S_{n} with dimension less than {dimension_limit}: {count}")

# To run this code, you might need to install the sympy library:
# pip install sympy
if __name__ == '__main__':
    solve_sn_irreps_dimension()
