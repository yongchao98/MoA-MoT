import math
from collections import Counter

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def solve_task():
    """
    Calculates the number of equivalence classes of endomorphisms on a set of size 4.
    This is equivalent to counting the number of non-isomorphic functional graphs on 4 vertices.
    The method is based on integer partitions and combinations with repetition.
    """
    n = 4
    
    # c[k] is the number of non-isomorphic connected functional graphs on k vertices.
    # These values are from OEIS A000248.
    c = {1: 1, 2: 2, 3: 4, 4: 9}

    # The integer partitions of 4 represent the sizes of the connected components.
    partitions_of_4 = [
        [4],
        [3, 1],
        [2, 2],
        [2, 1, 1],
        [1, 1, 1, 1]
    ]

    print(f"The number of equivalence classes of endomorphisms on a set of size {n} is found by summing the ways to form graphs for each integer partition of {n}:")
    print("-" * 75)
    
    results = []
    
    for p in partitions_of_4:
        counts = Counter(p)
        term = 1
        
        # The formula for a given partition is the product over distinct part sizes k of:
        # C(c[k] + m_k - 1, m_k) where m_k is the multiplicity of k in the partition.
        
        # Build a string for explaining the calculation
        calc_str_parts = []
        for k, m_k in sorted(counts.items(), reverse=True):
            c_k = c[k]
            term *= combinations(c_k + m_k - 1, m_k)
            if m_k == 1:
                calc_str_parts.append(f"c({k})")
            else:
                # This is the formula for combinations with repetition
                calc_str_parts.append(f"C(c({k})+{m_k}-1, {m_k})")
        
        calc_str = " * ".join(calc_str_parts)
        final_calc_str = calc_str
        for k_val in c:
            final_calc_str = final_calc_str.replace(f'c({k_val})', str(c[k_val]))

        print(f"For partition {str(p):<12}: {final_calc_str:<30} = {term}")
        results.append(term)

    total_classes = sum(results)
    
    print("-" * 75)
    # The final prompt wants each number in the final equation.
    equation_str = " + ".join(map(str, results))
    print(f"The total count is the sum of these values:")
    print(f"{equation_str} = {total_classes}")

solve_task()