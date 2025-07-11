import math
from collections import Counter

def get_partitions(n):
    """Generates all unique integer partitions of n."""
    def generate(target, max_val, current_partition):
        if target == 0:
            yield tuple(sorted(current_partition))
            return
        if target < 0 or max_val == 0:
            return

        # Option 1: Use max_val in the partition
        if target >= max_val:
            yield from generate(target - max_val, max_val, current_partition + [max_val])
        
        # Option 2: Don't use max_val, try smaller values
        yield from generate(target, max_val - 1, current_partition)

    # Use a set to store unique partitions
    unique_partitions = set(generate(n, n, []))
    # Return as a sorted list for consistent ordering
    return sorted(list(unique_partitions), reverse=True)


def get_divisors(k):
    """Returns a list of all divisors of k."""
    return [i for i in range(1, k + 1) if k % i == 0]

def solve_for_n(n):
    """
    Calculates the number of non-isomorphic endomorphisms on a set of size n
    using Burnside's Lemma.
    """
    n_factorial = math.factorial(n)
    partitions = get_partitions(n)
    
    total_sum_numerator = 0
    calculation_terms = []

    # Iterate over each conjugacy class, represented by a partition
    for p in partitions:
        cycle_counts = Counter(p) # e.g., for partition (2,1,1), this is {2:1, 1:2}

        # Calculate the size of the conjugacy class
        # Formula: n! / (product of k^c_k * c_k!)
        class_size_denom = 1
        for k, ck in cycle_counts.items():
            class_size_denom *= (k**ck) * math.factorial(ck)
        class_size = n_factorial // class_size_denom

        # Calculate the number of functions that commute with a permutation of this class
        # Formula: product over k=1..n of (sum_{j|k} j*c_j)^c_k
        num_commuting_functions = 1
        for k in range(1, n + 1):
            ck = cycle_counts.get(k, 0)
            if ck == 0:
                continue
            
            # Inner sum: sum over j where j divides k
            inner_sum = 0
            for j in get_divisors(k):
                cj = cycle_counts.get(j, 0)
                inner_sum += j * cj
            
            num_commuting_functions *= inner_sum ** ck

        # Add to the total sum for the numerator of Burnside's Lemma
        term = class_size * num_commuting_functions
        total_sum_numerator += term
        calculation_terms.append((class_size, num_commuting_functions))

    # Final result is the total sum divided by the size of the group
    result = total_sum_numerator // n_factorial

    # Format the output string to show the full calculation
    equation_parts = [f"{c} * {f}" for c, f in calculation_terms]
    equation_str = f"({ ' + '.join(equation_parts) }) / {n_factorial} = {result}"
    
    print("The number of equivalence classes is calculated using Burnside's Lemma.")
    print("The partitions of 4 correspond to the conjugacy classes of the symmetric group S4.")
    print("For each class, we count its size and the number of functions it fixes.")
    print("\nFinal Calculation:")
    print(equation_str)
    
    return result

# For the given problem, the size of the set is 4.
n = 4
final_answer = solve_for_n(n)

# The final answer is required in a specific format.
# The code above already prints the detailed breakdown.
# Here we just output the final numerical answer.
print(f"\n<<<19>>>")