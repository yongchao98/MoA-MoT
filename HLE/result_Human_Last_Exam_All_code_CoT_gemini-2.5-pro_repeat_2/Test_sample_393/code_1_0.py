import math
from collections import Counter

def solve_partition_count():
    """
    Calculates the number of ways to partition a set of n labeled items
    into unlabeled subsets of specified sizes.
    This corresponds to the number of Tverberg partitions under the problem's conditions.
    """
    # --- Example Input ---
    # n: total number of points
    # a: a list of integers representing the size of each part in the partition
    n = 8
    a = [3, 3, 2]

    # --- Verification (optional) ---
    if sum(a) != n:
        print(f"Error: The sum of parts {sum(a)} does not equal n {n}.")
        return

    # --- Step 1: Explain the formula ---
    print(f"The number of ways to partition a set of {n} items into subsets of sizes {a} is given by:")
    print(f"N = (1 / (product of factorials of counts of each size)) * (n! / (product of factorials of each part size))")
    print("-" * 20)

    # --- Step 2: Calculate components for the formula ---
    # Frequencies of each part size
    size_counts = Counter(a)
    m_values = list(size_counts.values())

    # --- Step 3: Display the formula with numbers ---
    n_str = f"{n}!"
    a_fact_str = " * ".join([f"{part}!" for part in a])
    m_fact_str = " * ".join([f"{count}!" for count in m_values])
    print("Plugging in the numbers for our example (n=8, parts=[3, 3, 2]):")
    print(f"N = (1 / ({m_fact_str})) * ({n_str} / ({a_fact_str}))")
    print("-" * 20)

    # --- Step 4: Calculate intermediate results ---
    n_factorial = math.factorial(n)
    a_factorials_prod = 1
    for part in a:
        a_factorials_prod *= math.factorial(part)

    m_factorials_prod = 1
    for count in m_values:
        m_factorials_prod *= math.factorial(count)

    print("Calculating the factorial values:")
    print(f"{n}! = {n_factorial}")
    a_fact_vals_str = " * ".join([str(math.factorial(part)) for part in a])
    print(f"Product of part size factorials = {a_fact_vals_str} = {a_factorials_prod}")
    m_fact_vals_str = " * ".join([str(math.factorial(count)) for count in m_values])
    print(f"Product of size count factorials = {m_fact_vals_str} = {m_factorials_prod}")
    print("-" * 20)

    # --- Step 5: Show the calculation steps ---
    print("Substituting these values into the formula:")
    print(f"N = (1 / {m_factorials_prod}) * ({n_factorial} / {a_factorials_prod})")

    multinomial_coeff = n_factorial // a_factorials_prod
    print(f"N = (1 / {m_factorials_prod}) * {multinomial_coeff}")

    final_result = multinomial_coeff // m_factorials_prod
    print(f"N = {final_result}")
    print("-" * 20)
    print(f"The final number of partitions is {final_result}.")


if __name__ == '__main__':
    solve_partition_count()
