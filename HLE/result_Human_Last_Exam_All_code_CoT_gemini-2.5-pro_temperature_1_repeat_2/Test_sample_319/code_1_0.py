import math

def get_partitions(n, min_part_size):
    """
    Generates all integer partitions of n where each part is at least min_part_size.
    """
    # Using a recursive generator
    if n == 0:
        yield []
        return
    
    # Start searching for parts from min_part_size up to n
    # The next part cannot be smaller than the previous one to ensure unique partitions
    def generate(target, start, current_partition):
        if target == 0:
            yield list(current_partition)
            return

        for i in range(start, target + 1):
            if i >= min_part_size:
                current_partition.append(i)
                yield from generate(target - i, i, current_partition)
                current_partition.pop()

    yield from generate(n, min_part_size, [])

def calculate_exponent(partition):
    """
    Calculates the sum of ceil(d_j / 2) for a given partition.
    """
    return sum(math.ceil(part / 2) for part in partition)

def solve():
    """
    Finds the largest possible value of c.
    """
    D = 10  # Dimension of the space
    m = 2   # Dimension of the planes' direction spaces
    min_part_size = m

    print(f"The problem is to find the maximum value of c.")
    print(f"This corresponds to maximizing the sum of ceil(d_j / 2) for all partitions d_j of {D},")
    print(f"where each part d_j must be at least {min_part_size}.")
    print("\nLet's analyze the possible partitions and their corresponding exponent c:")
    
    partitions = list(get_partitions(D, min_part_size))
    
    max_c = 0
    best_partitions = []
    
    for p in partitions:
        c = calculate_exponent(p)
        print(f"Partition: {str(p):<20} -> c = sum(ceil(d/2)) = {' + '.join([f'ceil({d}/2)' for d in p])} = {' + '.join([str(math.ceil(d/2)) for d in p])} = {c}")
        if c > max_c:
            max_c = c
            best_partitions = [p]
        elif c == max_c:
            best_partitions.append(p)
            
    print("\n-------------------------------------------------------------")
    print(f"The maximum possible value for c is {max_c}.")
    print(f"This maximum value is achieved for the following partitions:")
    for p in best_partitions:
        print(p)

solve()