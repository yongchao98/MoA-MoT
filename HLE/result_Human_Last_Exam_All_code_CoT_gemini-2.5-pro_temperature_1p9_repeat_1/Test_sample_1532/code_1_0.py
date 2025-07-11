import math

def solve_coalition_problem():
    """
    Calculates the smallest number of states in a "strong coalition"
    able to pass a constitutional amendment.
    """

    # The number of states in the USA is 50. D.C. and Puerto Rico do not
    # count as states for the purposes of Article V of the Constitution.
    total_states = 50

    # Step 1: Calculate the number of states needed to PROPOSE an amendment.
    # A coalition can bypass a hostile Congress by having 2/3 of state legislatures
    # call for a national convention.
    proposal_fraction = 2/3
    states_needed_to_propose = math.ceil(proposal_fraction * total_states)

    # Step 2: Calculate the number of states needed to RATIFY an amendment.
    # Ratification requires the approval of 3/4 of the states.
    ratification_fraction = 3/4
    states_needed_to_ratify = math.ceil(ratification_fraction * total_states)

    # Step 3: A strong coalition must meet the highest threshold to guarantee success.
    # The smallest size of such a coalition is the maximum of the numbers required
    # for proposal and ratification.
    smallest_strong_coalition = max(states_needed_to_propose, states_needed_to_ratify)

    # --- Output the step-by-step calculation ---
    print(f"Analysis of the Constitutional Amendment Process:")
    print(f"Total number of states = {total_states}")
    print("-" * 40)
    print("Step 1: Proposal Requirement")
    print(f"A proposal can be forced by 2/3 of state legislatures.")
    print(f"Calculation: ceil({proposal_fraction:.2f} * {total_states}) = ceil({proposal_fraction * total_states:.2f}) = {states_needed_to_propose} states")
    print("-" * 40)
    print("Step 2: Ratification Requirement")
    print(f"Ratification requires 3/4 of states.")
    print(f"Calculation: ceil({ratification_fraction:.2f} * {total_states}) = ceil({ratification_fraction * total_states:.2f}) = {states_needed_to_ratify} states")
    print("-" * 40)
    print("Step 3: Determining the Smallest Strong Coalition")
    print("A strong coalition must overcome the highest barrier (the step requiring the most states).")
    print(f"Final equation: Smallest coalition = max(states_to_propose, states_to_ratify)")
    print(f"Final numbers: Smallest coalition = max({states_needed_to_propose}, {states_needed_to_ratify})")
    print("-" * 40)
    print(f"The result is {smallest_strong_coalition}.")


solve_coalition_problem()
<<<38>>>