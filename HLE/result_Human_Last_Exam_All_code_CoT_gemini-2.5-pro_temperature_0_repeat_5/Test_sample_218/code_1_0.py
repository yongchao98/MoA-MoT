import math
from sympy.utilities.iterables import partitions

def calculate_dimension(partition):
    """
    Calculates the dimension of an irreducible representation of S_n
    corresponding to the given partition using the Hook Length Formula.
    """
    n = sum(partition)
    
    # Calculate the product of hook lengths
    hook_product = 1
    for i, part_len in enumerate(partition): # i is the row index (0-based)
        for j in range(part_len): # j is the column index (0-based)
            # Number of cells to the right
            hooks_right = part_len - (j + 1)
            
            # Number of cells below
            hooks_down = 0
            for k in range(i + 1, len(partition)):
                if partition[k] > j:
                    hooks_down += 1
            
            # Hook length is hooks_right + hooks_down + 1 (for the cell itself)
            hook_length = hooks_right + hooks_down + 1
            hook_product *= hook_length
            
    # Calculate dimension using the formula n! / product_of_hooks
    # Python's integers handle the large numbers involved
    if hook_product == 0:
        return 0
    dimension = math.factorial(n) // hook_product
    return dimension

def find_irreps_below_threshold():
    """
    Finds the number of irreps of S_25 with dimension less than a threshold.
    """
    n = 25
    threshold = 500000
    count = 0
    
    print(f"Finding irreps of S_{n} with dimension < {threshold}:\n")

    # Generate all partitions of n
    # The partitions() function from sympy yields dictionaries
    all_partitions_gen = partitions(n)

    for p_dict in all_partitions_gen:
        # Convert the partition from sympy's dict format to a tuple
        # e.g., {3:1, 1:2} becomes (3, 1, 1)
        p_list = []
        for val, num in sorted(p_dict.items(), reverse=True):
            p_list.extend([val] * num)
        
        partition_tuple = tuple(p_list)
        
        # Calculate the dimension for the current partition
        dim = calculate_dimension(partition_tuple)
        
        # Check if the dimension is below the threshold
        if dim < threshold:
            count += 1
            # As requested, output the numbers for each qualifying representation
            print(f"Partition {partition_tuple}: dimension = {dim}")

    print(f"\nTotal number of irreducible representations with dimension strictly less than {threshold}: {count}")

if __name__ == '__main__':
    find_irreps_below_threshold()
<<<62>>>