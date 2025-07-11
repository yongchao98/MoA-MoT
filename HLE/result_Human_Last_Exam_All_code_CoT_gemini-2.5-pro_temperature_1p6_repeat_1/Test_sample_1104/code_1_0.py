import math

def solve_proportionality_puzzle():
    """
    Calculates s1 and s2 based on the PJR and EJR conditions.
    """
    
    # Parameters from the problem
    k = 100  # Committee size
    
    # Analysis of the initial 6 voters
    # Group N' = {1, 2, 3, 4, 5, 6}
    l = 6  # Size of the group N'
    # Intersection of approved candidates for N' is {a, b, c}
    j = 3  # Size of the intersection
    
    # --- Part 1: Calculation of s1 (PJR) ---
    # The condition for a PJR-satisfying committee where voter 1 is unsatisfied is l < n/k.
    # We need to find the smallest integer n (which is s1) that satisfies this.
    # n > l * k
    s1 = l * k + 1
    
    print("--- Calculating s1 (for PJR) ---")
    print(f"A committee that leaves voter 1 unsatisfied must not contain candidates from {{a, b, c}}.")
    print(f"Consider the group of the first 6 voters. Its size is l = {l}.")
    print(f"Their commonly approved candidates are {{'a', 'b', 'c'}}, so j = {j}.")
    print("For PJR to be satisfied, we must avoid the case where a group is large enough (>= n/k) but has no representatives.")
    print(f"Since their common candidates are not in the committee, the group's size must be strictly less than the PJR threshold.")
    print(f"So, the inequality l < n/k must hold.")
    print(f"Substituting the values: {l} < n / {k}")
    print(f"This simplifies to: n > {l} * {k}, which is n > {l * k}.")
    print(f"The smallest integer n satisfying this is s1 = {s1}.\n")

    # --- Part 2: Calculation of s2 (EJR) ---
    # The condition for an EJR-satisfying committee where voter 1 is unsatisfied is l/n <= j/k.
    # We need to find the smallest integer n (which is s2) that satisfies this.
    # n >= l * k / j
    s2 = math.ceil(l * k / j)
    
    print("--- Calculating s2 (for EJR) ---")
    print(f"For EJR, we use the same group of {l} voters with {j} common candidates.")
    print("For EJR to be satisfied, the ratio of the group size to the total number of voters (l/n) must not be unjustly high compared to the ratio of their common candidates to the committee size (j/k).")
    print("Since their common candidates are not represented, the inequality l/n <= j/k must hold.")
    print(f"Substituting the values: {l} / n <= {j} / {k}")
    print(f"This simplifies to: n >= ({l} * {k}) / {j}")
    print(f"n >= {l * k} / {j} = {l * k / j}")
    print(f"The smallest integer n satisfying this is s2 = {s2}.\n")

    # --- Final Result ---
    result = (s1, s2)
    print(f"The solution pair (s1, s2) is: {result}")

solve_proportionality_puzzle()