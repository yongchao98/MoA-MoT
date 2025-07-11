import math

def solve_proportionality_riddles():
    """
    This script calculates the smallest preference profile sizes (s1 and s2)
    based on the definitions of PJR and EJR.
    """
    
    k = 100  # Committee size
    
    # --- Part 1: Calculation of s1 for PJR ---
    print("--- Calculating s1 (for PJR) ---")
    
    # The PJR rule states that for any cohesive group of unsatisfied voters N',
    # we must have |N'| < n/k.
    # We are given that voter 1 is unsatisfied. Voter 1 alone forms a cohesive
    # group of size 1. Let's call this group N_v1.
    size_N_v1 = 1
    
    # So, the inequality is |N_v1| < s1 / k
    print(f"The PJR condition for the group {{voter 1}} is: {size_N_v1} < s1 / {k}")
    
    # This implies s1 > size_N_v1 * k
    s1_lower_bound = size_N_v1 * k
    print(f"This simplifies to: s1 > {s1_lower_bound}")
    
    # The smallest integer s1 satisfying this is s1_lower_bound + 1
    s1 = s1_lower_bound + 1
    print(f"The smallest integer profile size s1 is therefore: {s1}")
    print("-" * 20)
    
    # --- Part 2: Calculation of s2 for EJR ---
    print("\n--- Calculating s2 (for EJR) ---")
    
    # The EJR rule states that for any group of unsatisfied voters N'
    # and any integer l >= 1 with |intersection(ballots)| >= l,
    # we must have |N'| <= l * n / k.
    # For the group N'={voter 1}, |N'|=1 and |A(1)|=4.
    # The condition must hold for l=1, 2, 3, and 4.
    
    size_N_v1 = 1
    ballot_size_A1 = 4
    
    print(f"The EJR condition for N'={{voter 1}} (size {size_N_v1}) and ballot intersection size {ballot_size_A1} must hold for l in [1, {ballot_size_A1}].")
    print("The condition is: |N'| <= l * s2 / k")

    # We check the condition for all possible l's to find the most constraining one.
    # s2 >= |N'| * k / l
    # The lower bound for s2 is max(|N'|*k/l for l=1..4). This is maximized when l is minimal (l=1).
    l_most_restrictive = 1
    
    print(f"\nThe most restrictive condition is for the smallest l, i.e., l = {l_most_restrictive}.")
    print(f"The inequality is: {size_N_v1} <= {l_most_restrictive} * s2 / {k}")
    
    # Calculate the lower bound for s2
    # s2 >= size_N_v1 * k / l_most_restrictive
    s2_lower_bound = (size_N_v1 * k) / l_most_restrictive
    print(f"This simplifies to s2 >= {s2_lower_bound}")
    
    # s2 is the smallest integer satisfying this.
    s2 = math.ceil(s2_lower_bound)
    print(f"The smallest integer profile size s2 is therefore: {s2}")
    print("-" * 20)

    # --- Final Result ---
    print(f"\nAs argued in the explanation, these lower bounds are achievable.")
    print(f"The final solution is the pair (s1, s2).")
    print(f"({s1}, {s2})")

solve_proportionality_riddles()