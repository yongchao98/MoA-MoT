import math

def solve_proportionality():
    """
    This script calculates the smallest preference profile sizes (s1, s2)
    based on the definitions of PJR and EJR and the problem's constraints.
    """
    
    # --- Problem Constants ---
    k = 100  # Committee size

    # --- Part 1: Calculation for s_1 (Proportional Justified Representation) ---
    
    # The critical group consists of voter 1 alone, so its size L is 1.
    # To satisfy PJR, the size of this unsatisfied group must be less than the threshold n/k.
    L1 = 1
    
    # The inequality is L < s_1 / k
    # s_1 > L1 * k
    s1 = math.floor(L1 * k) + 1

    # --- Part 2: Calculation for s_2 (Extended Justified Representation) ---

    # The critical group is voters {1, 2, 3, 4, 5, 6}, so its size L is 6.
    # To satisfy EJR, for j=1, the group's size must be less than the threshold n/k.
    L2 = 6
    
    # The inequality is L < s_2 / k
    # s_2 > L2 * k
    s2 = math.floor(L2 * k) + 1

    # --- Final Output ---

    print("Calculation for s_1 (PJR):")
    # Output each number in the final equation
    print(f"The critical inequality is L < n/k. With L={L1} and k={k}, we get:")
    print(f"{L1} < s_1 / {k}  =>  s_1 > {L1 * k}")
    print(f"The smallest integer s_1 is {s1}.")
    
    print("\nCalculation for s_2 (EJR):")
    # Output each number in the final equation
    print(f"The critical inequality is L < n/k (for j=1). With L={L2} and k={k}, we get:")
    print(f"{L2} < s_2 / {k}  =>  s_2 > {L2 * k}")
    print(f"The smallest integer s_2 is {s2}.")
    
    print(f"\nThus, the solution is the pair ({s1}, {s2}).")

solve_proportionality()