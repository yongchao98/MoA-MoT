import math

def solve_constitution_game():
    """
    Calculates the smallest number of states in a "strong coalition"
    capable of ratifying a constitutional amendment against opposition.
    """
    
    # Step 1: Define the number of states relevant to the amendment process.
    # Article V refers to "States". Washington D.C. and Puerto Rico are not states
    # for the purpose of this article.
    total_states = 50
    
    print(f"The analysis considers the {total_states} United States for the amendment process.")
    print("-" * 30)

    # Step 2: Calculate the number of states needed to propose an amendment
    # via a convention, bypassing Congress. This requires 2/3 of states.
    proposal_fraction_numerator = 2
    proposal_fraction_denominator = 3
    
    # We need the ceiling because you can't have a fraction of a state.
    # Even 33.1 states rounds up to 34.
    states_to_propose = math.ceil(total_states * proposal_fraction_numerator / proposal_fraction_denominator)
    
    print("To force a constitutional convention for proposing an amendment, a coalition needs 2/3 of the states.")
    print(f"Calculation: ceil({proposal_fraction_numerator}/{proposal_fraction_denominator} * {total_states}) = {states_to_propose}")
    print("-" * 30)

    # Step 3: Calculate the number of states needed to ratify a proposed amendment.
    # This requires 3/4 of states.
    ratification_fraction_numerator = 3
    ratification_fraction_denominator = 4

    # We need the ceiling here as well.
    states_to_ratify = math.ceil(total_states * ratification_fraction_numerator / ratification_fraction_denominator)
    
    print("To ratify a proposed amendment, a coalition needs 3/4 of the states.")
    print(f"Calculation: ceil({ratification_fraction_numerator}/{ratification_fraction_denominator} * {total_states}) = {states_to_ratify}")
    print("-" * 30)
    
    # Step 4: Determine the size of the "strong coalition".
    # The coalition must be large enough for the most difficult step, which is the one
    # requiring the most states.
    min_strong_coalition_size = max(states_to_propose, states_to_ratify)
    
    print("A 'strong coalition' must be able to both propose and ratify the amendment against opposition.")
    print("Therefore, the size of the coalition must be equal to the higher of the two thresholds.")
    print(f"Final Calculation: max({states_to_propose}, {states_to_ratify}) = {min_strong_coalition_size}")
    
if __name__ == "__main__":
    solve_constitution_game()
