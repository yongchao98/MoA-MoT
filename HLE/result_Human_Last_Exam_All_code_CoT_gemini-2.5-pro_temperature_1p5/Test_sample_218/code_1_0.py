import math

def generate_partitions(n, min_part=1):
    """
    Generates all integer partitions of n.
    Yields partitions with parts in non-decreasing order.
    """
    if n == 0:
        yield []
        return
    if n < min_part:
        return
    
    # Iterate through possible first parts `i`
    for i in range(min_part, n + 1):
        # Recursively find partitions for the remainder `n-i`
        for sub_partition in generate_partitions(n - i, i):
            yield [i] + sub_partition

def calculate_dimension_and_hooks(partition, n):
    """
    Calculates the dimension of an irrep of S_n using the Hook Length Formula.
    The input partition must be sorted in descending order.
    
    Returns a tuple containing:
    (dimension, numerator, denominator, list_of_hook_lengths)
    """
    numerator = math.factorial(n)
    denominator = 1
    hook_lengths = []
    
    for r in range(len(partition)):  # r is the 0-indexed row number
        for c in range(partition[r]):  # c is the 0-indexed column number
            
            # 1. Boxes to the right
            boxes_right = partition[r] - (c + 1)
            
            # 2. Boxes below
            boxes_below = 0
            for i in range(r + 1, len(partition)):
                if partition[i] > c:
                    boxes_below += 1
            
            # 3. The box itself (+1)
            hook_length = boxes_right + boxes_below + 1
            
            hook_lengths.append(hook_length)
            denominator *= hook_length

    # The dimension is guaranteed to be an integer
    dimension = numerator // denominator
    return dimension, numerator, denominator, hook_lengths

def solve():
    """
    Finds and counts the irreps of S_25 with dimension less than 500,000.
    """
    n = 25
    limit = 500000
    count = 0

    print(f"Finding irreducible representations of S_{n} with dimension < {limit}:\n")

    # Generate all partitions of n
    partition_generator = generate_partitions(n)
    
    for p_asc in partition_generator:
        # The hook length calculation expects parts in descending order
        partition = sorted(p_asc, reverse=True)
        
        dim, num, den, hooks = calculate_dimension_and_hooks(partition, n)
        
        if dim < limit:
            count += 1
            
            # Format the hook product string for readability
            hooks_str = "*".join(map(str, sorted(hooks, reverse=True)))
            
            # To avoid an extremely long line, we don't print the value of 25!
            # or the full hook product string if it's too long.
            if len(hooks_str) > 80:
                hooks_str_display = f"Product of {len(hooks)} hook lengths"
            else:
                hooks_str_display = hooks_str

            print(f"#{count}: Partition {str(partition):<45} -> Dimension = {dim}")
            print(f"     Calculation: {n}! / ({hooks_str_display})")
            print("-" * 60)


    print(f"\nTotal number of irreducible representations with dimension strictly less than {limit}: {count}")

if __name__ == '__main__':
    solve()
    print("\n<<<56>>>")