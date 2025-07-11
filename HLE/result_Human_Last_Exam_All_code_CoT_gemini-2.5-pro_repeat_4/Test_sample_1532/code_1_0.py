import math

def solve_constitution_game():
    """
    Calculates the smallest number of states in a "strong coalition" capable
    of passing a constitutional amendment against arbitrary opposition.
    """
    # According to Article V of the U.S. Constitution, amendments are proposed
    # and ratified based on the 50 states.
    total_states = 50

    # --- ANALYSIS ---
    # A "strong coalition" must be able to force an amendment through the entire
    # legal process despite opposition. We analyze the two-step process:
    # proposal and ratification.

    # STEP 1: PROPOSAL
    # To bypass a hostile Congress, the coalition can use the national convention
    # route, which requires 2/3 of state legislatures to make an application.
    proposal_fraction = 2/3
    states_for_proposal = math.ceil(proposal_fraction * total_states)

    # STEP 2: RATIFICATION
    # Once proposed, an amendment must be ratified by 3/4 of the states.
    ratification_fraction = 3/4
    states_for_ratification = math.ceil(ratification_fraction * total_states)

    # STEP 3: CONCLUSION
    # The coalition must be large enough to satisfy the stricter of the two
    # conditions. The size must be the maximum of the numbers required for
    # proposal and ratification.
    min_coalition_size = max(states_for_proposal, states_for_ratification)

    # --- OUTPUT ---
    print("Analyzing the U.S. Constitutional amendment process to find the size of a 'strong coalition'.")
    print("-" * 80)
    print(f"Total number of states considered for the amendment process: {total_states}")
    print("\nStep 1: States needed to PROPOSE an amendment via a national convention.")
    print(f"The requirement is 2/3 of the states.")
    print(f"Calculation: ceil(2/3 * {total_states}) = {states_for_proposal}")

    print("\nStep 2: States needed to RATIFY a proposed amendment.")
    print(f"The requirement is 3/4 of the states.")
    print(f"Calculation: ceil(3/4 * {total_states}) = {states_for_ratification}")

    print("\nStep 3: Determining the smallest strong coalition size.")
    print("The coalition must be able to perform both proposal and ratification on its own.")
    print("Therefore, its size must meet the larger of the two requirements.")
    print("\nFinal Equation:")
    print(f"Smallest Coalition Size = max(States to Propose, States to Ratify)")
    print(f"Smallest Coalition Size = max({states_for_proposal}, {states_for_ratification}) = {min_coalition_size}")

solve_constitution_game()