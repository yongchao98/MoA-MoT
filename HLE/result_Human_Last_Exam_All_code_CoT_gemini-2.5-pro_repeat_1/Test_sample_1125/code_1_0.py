import math

def solve_hat_puzzle():
    """
    Calculates the largest number of team members (N) who are guaranteed
    to determine their hat number using a multi-phase strategy.
    """

    # --- Initial Setup ---
    num_members = 12
    print(f"Solving for {num_members} team members.\n")

    # --- Phase 1: Four Groups of Three ---
    # The 12 members are split into 4 groups of 3. In each group, they form a
    # cycle of pairings. For a group of size k, ceil(k/2) members are guaranteed
    # to learn their number.
    group_size_1 = 3
    num_groups_1 = num_members // group_size_1
    knowns_per_group_1 = math.ceil(group_size_1 / 2.0)
    total_knowns_phase1 = num_groups_1 * knowns_per_group_1
    unknowns_after_phase1 = num_members - total_knowns_phase1

    print("Phase 1: Divide into 4 groups of 3 and use cyclic pairings.")
    print(f"Guaranteed 'knowns' per group: {int(knowns_per_group_1)}")
    print(f"Total 'knowns' from Phase 1: {int(total_knowns_phase1)}")
    print(f"Remaining 'unknowns': {int(unknowns_after_phase1)}\n")

    # --- Phase 2: One Group of Four ---
    # The 4 remaining 'unknowns' form a new group and repeat the strategy.
    group_size_2 = unknowns_after_phase1
    knowns_phase2 = math.ceil(group_size_2 / 2.0)
    unknowns_after_phase2 = group_size_2 - knowns_phase2

    print("Phase 2: The 4 'unknowns' form a new group with cyclic pairings.")
    print(f"Guaranteed new 'knowns' from Phase 2: {int(knowns_phase2)}")
    print(f"Remaining 'unknowns': {int(unknowns_after_phase2)}\n")

    # --- Phase 3: The Final Pair ---
    # The last 2 'unknowns' pair up. Since all other 10 numbers are known,
    # revealing one of their numbers allows the other to know theirs by elimination.
    knowns_phase3 = unknowns_after_phase2
    
    print("Phase 3: The final 2 'unknowns' pair up.")
    print("One's number is revealed, the other knows theirs by elimination.")
    print(f"Guaranteed new 'knowns' from Phase 3: {int(knowns_phase3)}\n")

    # --- Final Calculation ---
    total_guaranteed = total_knowns_phase1 + knowns_phase2 + knowns_phase3

    print("The total number of people guaranteed to know their number is the sum from all phases.")
    # The final equation is printed as requested
    print(f"Final Equation: {int(total_knowns_phase1)} + {int(knowns_phase2)} + {int(knowns_phase3)} = {int(total_guaranteed)}")
    
    print("\nThe largest possible value of N is:")
    print(int(total_guaranteed))


solve_hat_puzzle()
<<<12>>>