import math

def find_smallest_strong_coalition():
    """
    This script calculates the smallest number of states in a "strong coalition"
    capable of passing a constitutional amendment against arbitrary opposition.

    The logic follows the path of a national convention, which depends on the
    number of states, not their population.
    """

    # --- Introduction ---
    print("Step 1: Determine the easiest path for a coalition of a small number of states.")
    print("The path of proposal by a National Convention bypasses the population-based House of Representatives.")
    print("This path has two stages: proposal and ratification.\n")

    # --- Stage 1: Proposal by National Convention ---
    # The Constitution requires "two thirds of the several States" (the 50 states)
    # to call for a convention.
    num_states_for_proposal = 50
    proposal_threshold = math.ceil(2/3 * num_states_for_proposal)
    print(f"Step 2: Calculate the number of states needed for the proposal stage.")
    print(f"The coalition needs 2/3 of the {num_states_for_proposal} states to call a convention.")
    print(f"  ceil(2/3 * {num_states_for_proposal}) = {proposal_threshold} states.\n")

    # --- Stage 2: Ratification ---
    # Ratification requires "three fourths of the several States."
    # This includes the 50 states plus Washington D.C.
    num_entities_for_ratification = 51
    ratification_threshold = math.ceil(3/4 * num_entities_for_ratification)
    print(f"Step 3: Calculate the number of states/entities needed for ratification.")
    print(f"The coalition needs 3/4 of the {num_entities_for_ratification} ratifying entities (50 states + D.C.).")
    print("The final equation is based on this step, as it's the higher threshold:")
    # The final equation with each number, as requested.
    print(f"  3/4 * {num_entities_for_ratification} = {3/4 * num_entities_for_ratification}")
    print(f"  Taking the ceiling of that result gives the final number: {ratification_threshold}\n")

    # --- Conclusion ---
    strong_coalition_size = max(proposal_threshold, ratification_threshold)
    print("Step 4: Determine the minimum size of the strong coalition.")
    print(f"The coalition must meet the higher of the two thresholds ({proposal_threshold} and {ratification_threshold}).")
    print(f"The smallest number of states in a strong coalition is therefore {strong_coalition_size}.")


if __name__ == '__main__':
    find_smallest_strong_coalition()