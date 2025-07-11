import math

def solve_strong_coalition():
    """
    Calculates the smallest number of states in a "strong coalition"
    able to pass a constitutional amendment under adversarial conditions.
    """
    # 1. Identify the total number of "states" as per the problem description.
    total_states = 50 + 1  # 50 US States + Puerto Rico

    # 2. Define the constitutional fractions for proposal and ratification.
    proposal_fraction = 2/3
    ratification_fraction = 3/4

    # 3. Calculate the number of states needed for proposal via a convention.
    # This must be an integer, so we use ceiling to round up.
    states_for_proposal = math.ceil(proposal_fraction * total_states)

    # 4. Calculate the number of states needed for ratification.
    states_for_ratification = math.ceil(ratification_fraction * total_states)

    # 5. A "strong coalition" must meet the stricter (higher) of the two thresholds.
    smallest_coalition_size = max(states_for_proposal, states_for_ratification)

    # Print the step-by-step analysis.
    print(f"To find the smallest 'strong coalition', we analyze the constitutional amendment process.")
    print(f"First, we define our total number of states: 50 states + Puerto Rico = {total_states} states.")
    print("-" * 40)

    print("Step 1: The Proposal Stage")
    print("An amendment can be proposed by a convention called for by 2/3 of the states.")
    print(f"The number of states needed is ceil(2/3 * {total_states}) = {states_for_proposal}.")
    print("-" * 40)

    print("Step 2: The Ratification Stage")
    print("An amendment must be ratified by 3/4 of the states.")
    print(f"The number of states needed is ceil(3/4 * {total_states}) = {states_for_ratification}.")
    print("-" * 40)

    print("Conclusion: Finding the Smallest Strong Coalition")
    print(f"The coalition must have enough states to meet the requirements of both stages.")
    print(f"To guarantee success, they must satisfy the higher number.")
    print(f"The final calculation is: max({states_for_proposal}, {states_for_ratification}) = {smallest_coalition_size}")

if __name__ == '__main__':
    solve_strong_coalition()
    print("\n<<<39>>>")