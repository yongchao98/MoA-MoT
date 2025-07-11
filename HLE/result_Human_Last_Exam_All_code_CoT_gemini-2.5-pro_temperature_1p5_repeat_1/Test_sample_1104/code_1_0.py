import math

def solve_proportionality_problem():
    """
    This function calculates the smallest preference profile sizes (s1, s2)
    for PJR and EJR under the given conditions.
    """
    
    # --- Problem Parameters ---
    # k: The size of the committee.
    # V_size: The size of the critical group of voters {1, 2, 3, 4, 5, 6}.
    # intersection_size: The number of commonly approved candidates {a, b, c} by this group.
    k = 100
    V_size = 6
    intersection_size = 3

    # --- Part 1: Calculation for s1 (Proportional Justified Representation) ---
    print("--- Finding s1 for PJR ---")
    print("To leave voter 1 unsatisfied, candidates {a, b, c, x} cannot be in the committee W.")
    print(f"A critical group of voters is V = {{1, 2, 3, 4, 5, 6}}, with |V| = {V_size}.")
    print(f"This group is cohesive, with common approvals {{a, b, c}}, and is unsatisfied by W.")
    print("For PJR to hold, this group must not meet the size requirement for a justified protest.")
    print("The condition to avoid a PJR violation is: |V| < n / k")
    print(f"Substituting the values: {V_size} < n / {k}")
    pjr_threshold = V_size * k
    print(f"This implies: n > {V_size} * {k}, which is n > {pjr_threshold}.")
    s1 = math.floor(pjr_threshold) + 1
    print(f"The smallest integer n that satisfies this condition is {s1}.")
    print(f"Thus, s1 = {s1}\n")

    # --- Part 2: Calculation for s2 (Extended Justified Representation) ---
    print("--- Finding s2 for EJR ---")
    print(f"For EJR, we analyze the same critical group V, with |V| = {V_size}.")
    print(f"The size of the common intersection is |S| = {intersection_size}, so the group is l-cohesive for l=1, 2, 3.")
    print("To avoid an EJR violation, the condition |V| < l * n / k must hold for each l.")
    print("This can be rewritten as n > (|V| * k) / l.")

    max_ejr_threshold = 0
    for l in range(1, intersection_size + 1):
        ejr_threshold = (V_size * k) / l
        if ejr_threshold > max_ejr_threshold:
            max_ejr_threshold = ejr_threshold
        print(f"For l = {l}: n > ({V_size} * {k}) / {l}  =>  n > {ejr_threshold:.2f}")

    print(f"\nAll these conditions must be met. The strictest condition is n > {max_ejr_threshold:.2f}.")
    s2 = math.floor(max_ejr_threshold) + 1
    print(f"The smallest integer n that satisfies the strictest condition is {s2}.")
    print(f"Thus, s2 = {s2}\n")

    # --- Final Answer ---
    print(f"The final solution for the pair (s1, s2) is ({s1}, {s2}).")

# Execute the function to solve the problem
solve_proportionality_problem()