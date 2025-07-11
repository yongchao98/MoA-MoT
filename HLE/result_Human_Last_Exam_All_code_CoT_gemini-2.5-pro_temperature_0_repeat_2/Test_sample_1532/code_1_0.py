import math

def solve_strong_coalition():
    """
    Calculates the smallest number of states in a "strong coalition"
    capable of ratifying a constitutional amendment.
    """
    # The amendment process legally involves the 50 states.
    num_states = 50

    # Step 1: Calculate the number of states needed for proposal.
    # An amendment can be proposed by a national convention called for by 2/3 of the state legislatures.
    # This is the most direct route for a coalition of states to force a proposal against opposition.
    proposal_fraction = 2/3
    states_for_proposal = math.ceil(proposal_fraction * num_states)

    # Step 2: Calculate the number of states needed for ratification.
    # An amendment must be ratified by 3/4 of the states to become part of the Constitution.
    ratification_fraction = 3/4
    states_for_ratification = math.ceil(ratification_fraction * num_states)

    # Step 3: Determine the size of the "strong coalition".
    # A strong coalition must be able to succeed at both proposing and ratifying the amendment.
    # Therefore, it must be large enough to meet the more stringent requirement (the maximum of the two).
    smallest_coalition_size = max(states_for_proposal, states_for_ratification)

    # Print the step-by-step reasoning and the final equation.
    print("Analyzing the Constitutional Amendment process (Article V):")
    print(f"Total number of states considered for amendment ratification: {num_states}")
    print("-" * 60)
    
    print("Step 1: Proposal Requirement")
    print("An amendment can be proposed by a convention called for by 2/3 of the states.")
    print(f"Calculation: ceil(2/3 * {num_states}) = {states_for_proposal}")
    print(f"Result: {states_for_proposal} states are needed to force a proposal.")
    print("-" * 60)

    print("Step 2: Ratification Requirement")
    print("An amendment must be ratified by 3/4 of the states.")
    print(f"Calculation: ceil(3/4 * {num_states}) = {states_for_ratification}")
    print(f"Result: {states_for_ratification} states are needed to ensure ratification.")
    print("-" * 60)

    print("Step 3: Finding the Smallest Strong Coalition")
    print("A 'strong coalition' must succeed despite opposition, so it must be large")
    print("enough to control the most demanding step of the process.")
    print("\nThe final equation is to find the maximum of the two requirements:")
    print(f"Smallest number of states = max({states_for_proposal}, {states_for_ratification}) = {smallest_coalition_size}")

solve_strong_coalition()