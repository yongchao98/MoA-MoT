import math
from collections import Counter

def n_choose_k(n, k):
    """Helper function for combinations"""
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def get_partitions(n):
    """Generates all integer partitions of n"""
    a = [0] * (n + 1)
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield tuple(a[:k + 2])
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield tuple(a[:k + 1])

def main():
    """
    Calculates the number of equivalence classes of endomorphisms on a set of size 4.
    This is equivalent to counting non-isomorphic functional digraphs on 4 vertices.
    """
    n = 4
    
    # D_c(k) is the number of connected non-isomorphic functional digraphs on k vertices.
    # These values are known from combinatorial literature (e.g., Harary, 1959 or OEIS A000248).
    Dc = {1: 1, 2: 2, 3: 4, 4: 10}
    
    total_classes = 0
    
    print(f"Calculating the number of equivalence classes for n = {n}\n")
    print(f"Using pre-computed counts of connected functional digraphs (Dc):")
    for k, v in Dc.items():
        print(f"Dc({k}) = {v}")
    print("-" * 30)

    partitions = list(get_partitions(n))
    
    for p in partitions:
        counts = Counter(p)
        term_result = 1
        
        # Calculate combinations with repetition for each part size in the partition
        for size, multiplicity in counts.items():
            num_types = Dc.get(size, 0)
            term_result *= n_choose_k(num_types + multiplicity - 1, multiplicity)
            
        print(f"Partition {p}: ", end="")
        
        # Explain the calculation for the current partition
        parts_str = []
        for size, multiplicity in sorted(counts.items()):
            if multiplicity == 1:
                parts_str.append(f"Dc({size})")
            else:
                parts_str.append(f"multiset of size {multiplicity} from {Dc[size]} types")
        print(" * ".join(parts_str), end="")
        
        # Show the formula with numbers
        formula_parts = []
        for size, multiplicity in sorted(counts.items()):
            num_types = Dc.get(size, 0)
            formula_parts.append(f"C({num_types}+{multiplicity-1}, {multiplicity})")
        print(f" = {' * '.join(formula_parts)} = {term_result}")

        total_classes += term_result

    print("-" * 30)
    print(f"The total number of classes is the sum of the results for all partitions.")
    partition_sums = []
    # Re-calculate to show the final equation
    for p in partitions:
        counts = Counter(p)
        term_result = 1
        for size, multiplicity in counts.items():
            num_types = Dc.get(size, 0)
            term_result *= n_choose_k(num_types + multiplicity - 1, multiplicity)
        partition_sums.append(str(term_result))
    
    final_equation = " + ".join(partition_sums)
    print(f"Final Sum: {final_equation} = {total_classes}")

if __name__ == "__main__":
    main()
