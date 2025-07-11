import math
from collections import Counter

def calculate_tverberg_partitions():
    """
    Calculates the number of ways to partition a set of n points into subsets
    of specified sizes, which corresponds to the number of Tverberg partitions
    for the specific set X described in the problem.
    """
    # --- User Input ---
    # Define the partition a_1, a_2, ..., a_r here.
    # This is a list of integers representing the sizes of the parts.
    # Example: For a partition of n=7 into parts of size 3, 2, 2, use [3, 2, 2].
    partition_parts = [3, 2, 2]

    # --- Calculation ---

    # 1. Calculate n, the total number of points.
    n = sum(partition_parts)
    print(f"The given partition is: {partition_parts}")
    print(f"The total number of points, n, is the sum of the parts:")
    # Create the string for the sum equation
    sum_equation = f"n = {' + '.join(map(str, partition_parts))} = {n}"
    print(sum_equation)
    print("-" * 30)

    # 2. Calculate the multinomial coefficient for distinguishable parts.
    n_factorial = math.factorial(n)
    a_factorials = [math.factorial(p) for p in partition_parts]
    denominator_a_factorial = math.prod(a_factorials)
    multinomial_coefficient = n_factorial // denominator_a_factorial

    print("Step 1: Calculate the number of ways to create ordered partitions (multinomial coefficient).")
    print(f"Formula: n! / (a_1! * a_2! * ...)")
    # Create strings for the equation steps
    a_factorial_str = ' * '.join([f"{p}!" for p in partition_parts])
    a_factorial_values_str = ' * '.join(map(str, a_factorials))
    print(f"Equation: {n}! / ({a_factorial_str})")
    print(f"= {n_factorial} / ({a_factorial_values_str})")
    print(f"= {n_factorial} / {denominator_a_factorial}")
    print(f"= {multinomial_coefficient}")
    print("-" * 30)

    # 3. Find counts of identical part sizes (k_j) to correct for indistinguishable parts.
    part_size_counts = Counter(partition_parts)
    counts_k = list(part_size_counts.values())
    k_factorials = [math.factorial(c) for c in counts_k]
    denominator_k_factorial = math.prod(k_factorials)

    print("Step 2: Correct for indistinguishable parts (parts of the same size).")
    print("We count how many times each part size appears:")
    for size, count in part_size_counts.items():
        print(f"  - Part size {size} appears {count} time(s).")
    
    print("We must divide by the factorial of these counts.")
    # Create strings for the correction factor equation
    k_factorial_str = ' * '.join([f"{c}!" for c in counts_k])
    k_factorial_values_str = ' * '.join(map(str, k_factorials))
    print(f"Correction Factor = {k_factorial_str}")
    print(f"= {k_factorial_values_str}")
    print(f"= {denominator_k_factorial}")
    print("-" * 30)

    # 4. Calculate the final result.
    final_result = multinomial_coefficient // denominator_k_factorial

    print("Step 3: Calculate the final number of Tverberg partitions.")
    print("Result = (Multinomial Coefficient) / (Correction Factor)")
    print(f"Result = {multinomial_coefficient} / {denominator_k_factorial}")
    print(f"Result = {final_result}")
    print("=" * 30)

    # Final summary equation
    print("The final calculation can be expressed as:")
    final_equation_str = f"({n}! / ({a_factorial_str})) / ({k_factorial_str}) = {final_result}"
    print(final_equation_str)

if __name__ == '__main__':
    calculate_tverberg_partitions()
<<<105>>>