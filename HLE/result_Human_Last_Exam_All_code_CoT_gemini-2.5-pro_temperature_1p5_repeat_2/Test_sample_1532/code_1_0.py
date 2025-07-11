import math

def solve_strong_coalition():
    """
    Calculates the smallest number of states in a "strong coalition"
    capable of ratifying a constitutional amendment against arbitrary opposition.
    """

    num_states = 50
    proposal_fraction = 2/3
    ratification_fraction = 3/4

    print("Step 1: Determine the number of states needed for the proposal stage.")
    print("In an adversarial scenario, the coalition must bypass Congress.")
    print("This requires 2/3 of state legislatures to call a convention.")
    
    # Calculate the number of states needed for proposal
    proposal_needed = math.ceil(num_states * proposal_fraction)
    print(f"Calculation: ceil({proposal_fraction:.2f} * {num_states}) = ceil({num_states * proposal_fraction:.2f}) = {proposal_needed}\n")

    print("Step 2: Determine the number of states needed for the ratification stage.")
    print("An amendment must be ratified by 3/4 of the states.")

    # Calculate the number of states needed for ratification
    ratification_needed = math.ceil(num_states * ratification_fraction)
    print(f"Calculation: ceil({ratification_fraction:.2f} * {num_states}) = ceil({num_states * ratification_fraction:.2f}) = {ratification_needed}\n")

    print("Step 3: Determine the size of the 'strong coalition'.")
    print("The coalition must be large enough to clear the higher of the two hurdles (proposal and ratification).")
    
    # The size of the strong coalition is the maximum of the two requirements
    strong_coalition_size = max(proposal_needed, ratification_needed)
    
    print(f"Final Equation: Smallest Coalition Size = max({proposal_needed}, {ratification_needed})")
    print(f"The smallest number of states that could form a strong coalition is {strong_coalition_size}.")


solve_strong_coalition()
<<<38>>>