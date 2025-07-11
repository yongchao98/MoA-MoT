import math
from collections import Counter

def get_partitions(n):
    """
    Generates all integer partitions of n.
    """
    if n == 0:
        yield []
        return
    for p in get_partitions(n - 1):
        yield [1] + p
        if p and (len(p) < 2 or p[1] > p[0]):
            yield [p[0] + 1] + p[1:]

def combinations(n, k):
    """
    Calculates the number of combinations "n choose k".
    """
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve():
    """
    Calculates the number of non-isomorphic endomorphisms on a set of size 4.
    """
    n = 4

    # c(k) is the number of non-isomorphic CONNECTED functional graphs on k vertices.
    # From OEIS A001372: c(1)=1, c(2)=2, c(3)=4, c(4)=9.
    c = {1: 1, 2: 2, 3: 4, 4: 9}

    partitions = sorted(list(get_partitions(n)))
    
    total_count = 0
    contributions = []

    print(f"Calculating the number of non-isomorphic endomorphisms for a set of size {n}.")
    print("This is done by summing contributions from each integer partition of 4.\n")
    
    for p in partitions:
        # Count the frequency of each part size in the partition
        freq = Counter(p)
        
        contribution = 1
        calc_str_parts = []
        for k, j in freq.items():
            # For each part size k appearing j times, we choose j connected graphs of size k
            # from c(k) available types, with replacement.
            # The number of ways is C(c(k) + j - 1, j).
            term = combinations(c[k] + j - 1, j)
            contribution *= term
            calc_str_parts.append(f"C({c[k]}+{j}-1,{j}) for k={k}")

        contributions.append(contribution)
        total_count += contribution
        
        print(f"Partition {p}:")
        print(f"  Calculation: {' * '.join(calc_str_parts)}")
        print(f"  Contribution = {contribution}\n")

    # Format the final equation string
    equation_str = " + ".join(map(str, contributions))
    
    print("The total number is the sum of these contributions:")
    print(f"{equation_str} = {total_count}")

solve()
<<<19>>>