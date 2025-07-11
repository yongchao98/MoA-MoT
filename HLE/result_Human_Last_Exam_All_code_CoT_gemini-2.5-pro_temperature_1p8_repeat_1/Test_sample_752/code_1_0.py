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
    cycle structure (partition).
    """
    counts = Counter(partition)
    
    denominator = 1
    for k, m_k in counts.items():
        denominator *= (k**m_k) * math.factorial(m_k)
    
    # Size of conjugacy class = n! / denominator
    # Size of centralizer = n! / class_size = denominator
    return denominator

def solve_for_n(n):
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
    equation_str = " + ".join(map(str, sorted(centralizer_sizes)))
    print(f"The total cardinality is the sum of centralizer sizes for each of the {len(unique_partitions)} conjugacy classes of S_{n}:")
    print(f"\nEquation: {equation_str}")
    print(f"\nTotal Cardinality = {total_cardinality}")
    
    return total_cardinality

if __name__ == '__main__':
    n = 9
    result = solve_for_n(n)
    # The final answer format is specified by the problem.
    # The final print is not needed if the function solve_for_n already prints it.
