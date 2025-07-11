import math
from collections import Counter

def calculate_partitions():
    """
    Calculates the number of ways to partition a set of n points into
    unlabeled subsets of specified sizes.

    The problem states that for a given partition of n, a_1, ..., a_r, there
    exists a set of n points X for which *every* Tverberg partition has part
    sizes corresponding to a_1, ..., a_r. This simplifies the problem to a
    combinatorial one: counting the number of ways to partition the n labeled
    points of X into unlabeled subsets with the given sizes.

    The formula is: (n! / (a_1! * a_2! * ...)) / (m_1! * m_2! * ...),
    where n is the total number of points, a_i are the sizes of the parts,
    and m_j is the count of parts of a certain size.
    """

    # Example: A partition of n points into parts of given sizes.
    # For this to be a valid Tverberg partition in R^d, each part size a_i
    # must be at most d+1. For a = [4, 2, 2, 1], we would need d >= 3.
    a = [4, 2, 2, 1]
    
    n = sum(a)
    r = len(a)

    print(f"Given the partition of n into r={r} parts with sizes: {a}")
    print(f"The total number of points is n = sum({a}) = {n}\n")

    # --- Step 1: Calculate the multinomial coefficient ---
    # This counts partitions into LABELED/DISTINCT groups of sizes a_i.
    # Formula: n! / (a_1! * a_2! * ... * a_r!)

    n_factorial = math.factorial(n)
    
    # Calculate the product of factorials of each part size a_i
    denom_a_factorials_vals = [math.factorial(i) for i in a]
    denom_a_factorials_prod = math.prod(denom_a_factorials_vals)
    
    multinomial_coeff = n_factorial // denom_a_factorials_prod

    # Print the formula and calculation for the multinomial coefficient
    a_fact_str = " * ".join([f"{x}!" for x in a])
    denom_a_fact_vals_str = " * ".join(map(str, denom_a_factorials_vals))
    print("Step 1: Calculate ways to partition n items into LABELED groups of specified sizes.")
    print(f"Formula: n! / ({a_fact_str})")
    print(f"Calculation: {n}! / ({denom_a_fact_vals_str}) = {n_factorial} / {denom_a_factorials_prod} = {multinomial_coeff}\n")

    # --- Step 2: Correct for UNLABELED/INDISTINGUISHABLE groups ---
    # We divide by the permutations of groups of the same size.
    # This is equivalent to m_1! * m_2! * ..., where m_j is the frequency
    # of each distinct group size.

    size_counts = Counter(a)
    print(f"The sizes of the parts are {list(size_counts.keys())}.")
    print(f"The number of times each size appears is {list(size_counts.values())}.")

    m_factorials_vals = [math.factorial(count) for count in size_counts.values()]
    divisor_for_unlabeled = math.prod(m_factorials_vals)
    
    # Print the formula and calculation for the correction factor
    m_fact_str = " * ".join([f"{count}!" for count in size_counts.values()])
    m_fact_vals_str = " * ".join(map(str, m_factorials_vals))
    print("\nStep 2: Correct for UNLABELED groups by dividing by the product of factorials of size counts.")
    print(f"Divisor = {m_fact_str}")
    print(f"Calculation: {m_fact_vals_str} = {divisor_for_unlabeled}\n")

    # --- Step 3: Final Calculation ---
    final_result = multinomial_coeff // divisor_for_unlabeled

    print("Step 3: Final result is the value from Step 1 divided by the value from Step 2.")
    print(f"Total Number of Tverberg Partitions = {multinomial_coeff} / {divisor_for_unlabeled} = {final_result}")
    
    # Return the final number for the answer block
    return final_result

if __name__ == '__main__':
    result = calculate_partitions()
    print(f"\n<<< {result} >>>")
