import math

def solve_strong_coalition():
    """
    Calculates the smallest number of states in a "strong coalition"
    capable of passing a constitutional amendment.
    """

    # As per Article V of the US Constitution, only the 50 states can
    # participate in the amendment process.
    num_states = 50

    # --- Proposal Requirement ---
    # The coalition must be able to force a national convention for proposing
    # amendments. This requires an application from 2/3 of state legislatures.
    proposal_fraction = 2/3
    states_for_proposal = math.ceil(proposal_fraction * num_states)

    # --- Ratification Requirement ---
    # The coalition must be able to ratify the amendment. This requires
    # approval from 3/4 of the states.
    ratification_fraction = 3/4
    states_for_ratification = math.ceil(ratification_fraction * num_states)

    # A "strong coalition" must be able to clear the highest hurdle on its own.
    # Therefore, its size must be the maximum of the two requirements.
    min_coalition_size = max(states_for_proposal, states_for_ratification)

    print("Step 1: Determine the number of states needed to propose an amendment via convention.")
    print(f"  Calculation: ceil(2/3 * {num_states}) = {states_for_proposal}")
    print("-" * 30)

    print("Step 2: Determine the number of states needed to ratify an amendment.")
    print(f"  Calculation: ceil(3/4 * {num_states}) = {states_for_ratification}")
    print("-" * 30)

    print("Step 3: Determine the minimum size of a 'strong coalition'.")
    print("A strong coalition must be large enough to meet the higher of the two thresholds.")
    print(f"  Final Equation: Smallest Coalition Size = max({states_for_proposal}, {states_for_ratification})")
    print(f"\nThe smallest number of the mentioned States that could form a strong coalition is {min_coalition_size}.")


solve_strong_coalition()