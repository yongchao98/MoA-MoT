import math
from collections import Counter

def solve_partition_problem():
    """
    Calculates the number of Tverberg partitions based on the problem's premises.
    The problem is interpreted as counting the number of set partitions with specific part sizes.
    """

    # We select a concrete example for the partition a = [a_1, ..., a_r].
    # As per the problem description, these should satisfy certain constraints related
    # to Tverberg's theorem, but the calculation itself is purely combinatorial.
    # Let's use the partition a = [3, 3, 2, 2].
    # This implies n=10 elements partitioned into r=4 groups.
    # For context, if d=2, this partition is valid as n=10 >= (4-1)(2+1)+1=10
    # and all a_i <= d+1=3.
    a = [3, 3, 2, 2]
    n = sum(a)

    # --- Start of printing the explanation and calculation ---
    
    print("The number of Tverberg partitions for the given conditions is equivalent to counting")
    print(f"the number of ways to partition a set of {n} elements into unlabeled subsets of sizes {a}.")
    print("-" * 40)
    print("The general formula is: n! / ( (Π a_i!) * (Π k_j!) )")
    print("where a_i are the part sizes and k_j are the counts of identical part sizes.")
    print("-" * 40)

    # Calculate all necessary components for the formula
    
    # 1. Numerator: n!
    n_factorial = math.factorial(n)

    # 2. Denominator part 1: product of a_i!
    prod_a_factorial = 1
    for part_size in a:
        prod_a_factorial *= math.factorial(part_size)

    # 3. Denominator part 2: product of k_j! (from counts of identical part sizes)
    counts = Counter(a)
    prod_counts_factorial = 1
    for count in counts.values():
        prod_counts_factorial *= math.factorial(count)

    # 4. Final Result
    result = n_factorial // (prod_a_factorial * prod_counts_factorial)
    
    # --- Print the full equation with numbers ---
    
    # Build strings for the equation
    a_fact_str = " * ".join([f"{x}!" for x in a])
    a_fact_vals_str = " * ".join([str(math.factorial(x)) for x in a])
    
    counts_list = list(counts.values())
    counts_fact_str = " * ".join([f"{x}!" for x in counts_list])
    counts_fact_vals_str = " * ".join([str(math.factorial(x)) for x in counts_list])

    print("Calculation for the specific partition a = [3, 3, 2, 2]:")
    print(f"n = {n}")
    
    # Final equation format requested by the user
    print("\n--- Final Equation ---")
    equation_str_1 = f"Number of Partitions = {n}! / ( ({a_fact_str}) * ({counts_fact_str}) )"
    print(equation_str_1)
    
    equation_str_2 = f"                     = {n_factorial} / ( ({a_fact_vals_str}) * ({counts_fact_vals_str}) )"
    print(equation_str_2)
    
    equation_str_3 = f"                     = {n_factorial} / ( {prod_a_factorial} * {prod_counts_factorial} )"
    print(equation_str_3)

    equation_str_4 = f"                     = {n_factorial} / {prod_a_factorial * prod_counts_factorial}"
    print(equation_str_4)
    
    equation_str_5 = f"                     = {result}"
    print(equation_str_5)
    
    return result

# Execute the function
final_answer = solve_partition_problem()

# The final answer in the required format
print(f"<<<{final_answer}>>>")