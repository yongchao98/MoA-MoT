import math
from collections import Counter

def partitions_recursive(n, max_val=None):
    """
    Generates all integer partitions of n.
    Each partition is a list of integers that sum to n.
    """
    if max_val is None:
        max_val = n
    if n == 0:
        yield []
        return
    for i in range(min(n, max_val), 0, -1):
        for p in partitions_recursive(n - i, i):
            yield [i] + p

def get_centralizer_size(n, partition):
    """
    Calculates the size of the centralizer for an element with a given
    cycle structure (partition). This size is given by the formula
    z(λ) = product_k k^m_k * m_k! for a partition λ = (1^m_1 2^m_2 ...).
    """
    counts = Counter(partition)
    
    centralizer_size = 1
    for k, m_k in counts.items():
        centralizer_size *= (k**m_k) * math.factorial(m_k)
    
    return centralizer_size

def solve(n):
    """
    Calculates the cardinality of Γ(F)(•,•) for Σ_n.
    This is the sum of the sizes of centralizers for each conjugacy class.
    """
    unique_partitions = list(partitions_recursive(n))
    
    centralizer_sizes = []
    for part in unique_partitions:
        size = get_centralizer_size(n, part)
        centralizer_sizes.append(size)
        
    total_cardinality = sum(centralizer_sizes)
    
    # Output the full equation as requested
    equation_str_parts = [str(s) for s in sorted(centralizer_sizes)]
    print(f"The desired cardinality is the sum of the centralizer sizes for each of the {len(unique_partitions)} conjugacy classes of the symmetric group S_{n}.")
    print(f"For n={n}, the centralizer sizes are:")
    # print each number
    for size in sorted(centralizer_sizes):
        print(size)

    print("\nThe sum is:")
    print(" + ".join(equation_str_parts), "=", total_cardinality)
    
    return total_cardinality

if __name__ == '__main__':
    solve(9)
