import math

def solve_strong_coalition():
    """
    Calculates the smallest number of states in a "strong coalition"
    capable of ratifying a constitutional amendment.
    """
    # Step 1: Define the total number of entities that can participate in the amendment process.
    # According to Article V of the U.S. Constitution, only "States" can ratify amendments.
    # The District of Columbia and Puerto Rico are not states for this purpose.
    total_states = 50
    print(f"The analysis is based on the {total_states} U.S. States, as only they can participate in the amendment process under Article V.")
    print("-" * 30)

    # Step 2: Calculate the number of states needed to propose an amendment.
    # The most direct route for a coalition of states is to have 2/3 of state legislatures
    # call for a national convention, bypassing the U.S. Congress.
    proposal_fraction = 2/3
    states_for_proposal = math.ceil(proposal_fraction * total_states)
    print("To guarantee a proposal, a coalition needs 2/3 of the states to call for a convention.")
    print(f"Number of states needed to propose: ceil(2/3 * {total_states}) = {states_for_proposal}")
    print("-" * 30)

    # Step 3: Calculate the number of states needed to ratify an amendment.
    # Ratification requires approval from 3/4 of the states.
    ratification_fraction = 3/4
    states_for_ratification = math.ceil(ratification_fraction * total_states)
    print("To guarantee ratification, the amendment must be approved by 3/4 of the states.")
    print(f"Number of states needed to ratify: ceil(3/4 * {total_states}) = {states_for_ratification}")
    print("-" * 30)

    # Step 4: Determine the size of the "strong coalition".
    # A strong coalition must be able to both propose and ratify. Therefore, it must
    # have enough states to clear the higher of the two thresholds.
    strong_coalition_size = max(states_for_proposal, states_for_ratification)
    print("A 'strong coalition' must be large enough to meet the ratification requirement, as it is the higher threshold.")
    print("\nFinal Equation:")
    print(f"max({states_for_proposal}, {states_for_ratification}) = {strong_coalition_size}")

solve_strong_coalition()