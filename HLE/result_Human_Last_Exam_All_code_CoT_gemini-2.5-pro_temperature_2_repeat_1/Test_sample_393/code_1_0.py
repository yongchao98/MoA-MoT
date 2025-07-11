import math
from collections import Counter

def solve_tverberg_partitions():
    """
    Calculates the number of Tverberg partitions for a given partition of n.

    The problem statement implies that the number of Tverberg partitions for a specially
    constructed point set X is equal to the number of ways to partition a set of n elements
    into subsets with specified sizes. This function calculates this combinatorial quantity.
    """
    # As per the plan, we use the example partition [3, 3, 2, 2].
    parts = [3, 3, 2, 2]
    n = sum(parts)
    
    # Numerator of the multinomial coefficient: n!
    numerator = math.factorial(n)
    
    # Denominator part 1: product of a_i! for each part a_i
    part_factorials_prod = 1
    part_factorials_str = []
    for p in parts:
        f = math.factorial(p)
        part_factorials_prod *= f
        part_factorials_str.append(f"{p}!")

    # Denominator part 2: correction for identical part sizes
    # Count frequencies of each part size
    counts = Counter(parts)
    
    counts_factorial_prod = 1
    counts_factorial_str = []
    for c in counts.values():
        f = math.factorial(c)
        counts_factorial_prod *= f
        counts_factorial_str.append(f"{c}!")

    # Calculate final result
    result = numerator // (part_factorials_prod * counts_factorial_prod)
    
    # Print the explanation and the step-by-step calculation
    print(f"The number of ways to partition a set of {n} items into parts of sizes {parts} is calculated as:")
    print(f"Formula: n! / (a1! * a2! * ...) / (c1! * c2! * ...)")
    print()
    print("Step 1: Calculate the multinomial coefficient for ordered partitions.")
    print(f"  {n}! / ({' * '.join(part_factorials_str)})")
    print(f"= {numerator} / ({' * '.join(map(str, [math.factorial(p) for p in parts]))})")
    multinomial_coeff = numerator // part_factorials_prod
    print(f"= {numerator} / {part_factorials_prod}")
    print(f"= {multinomial_coeff}")
    print()
    print("Step 2: Correct for unordered partitions by dividing by factorials of size counts.")
    print(f"The part sizes {parts} have repetitions. The counts of distinct sizes are {dict(counts)}.")
    print(f"Correction factor: {' * '.join(counts_factorial_str)} = {' * '.join(map(str, [math.factorial(c) for c in counts.values()]))} = {counts_factorial_prod}")
    print()
    print("Final Calculation:")
    print(f"{multinomial_coeff} / {counts_factorial_prod} = {result}")
    
solve_tverberg_partitions()
<<<6300>>>