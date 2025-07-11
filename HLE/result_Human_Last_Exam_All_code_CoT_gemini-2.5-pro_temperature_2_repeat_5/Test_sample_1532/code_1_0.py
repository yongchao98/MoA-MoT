import math

# Plan:
# 1. Define the number of states that can ratify amendments (50).
# 2. Calculate the number of states needed to meet the proposal requirement (2/3 of states).
# 3. Calculate the number of states needed to meet the ratification requirement (3/4 of states).
# 4. The size of a "strong coalition" is the maximum of these two numbers, as it must be able to
#    succeed at both stages.

# The number of states in the USA that can participate in the amendment process as of the year 2000.
# While DC and Puerto Rico are mentioned, they are not "States" in the sense of Article V.
num_states = 50
print(f"The analysis is based on the {num_states} states eligible to amend the Constitution.")
print("-" * 60)

# --- Stage 1: Proposal Requirement ---
print("Analyzing the Proposal Stage:")
print("A coalition of states can force an amendment proposal by having 2/3 of state legislatures call for a national convention.")
print("This path depends only on the number of states in the coalition.")

proposal_numerator = 2
proposal_denominator = 3
# We need to round up to the nearest whole number of states, so we use math.ceil()
states_for_proposal = math.ceil(num_states * proposal_numerator / proposal_denominator)

# Output the equation and result
print(f"Minimum states to propose = ceil(({proposal_numerator}/{proposal_denominator}) * {num_states}) = ceil({(num_states * proposal_numerator / proposal_denominator):.2f}) = {states_for_proposal}")
print("-" * 60)


# --- Stage 2: Ratification Requirement ---
print("Analyzing the Ratification Stage:")
print("After proposal, an amendment must be ratified by 3/4 of the states to become law.")

ratification_numerator = 3
ratification_denominator = 4
states_for_ratification = math.ceil(num_states * ratification_numerator / ratification_denominator)

# Output the equation and result
print(f"Minimum states to ratify = ceil(({ratification_numerator}/{ratification_denominator}) * {num_states}) = ceil({(num_states * ratification_numerator / ratification_denominator):.2f}) = {states_for_ratification}")
print("-" * 60)

# --- Conclusion ---
print("Determining the Strong Coalition Size:")
print("A 'strong coalition' must have enough states to overcome the hurdles of both proposal AND ratification.")
print("Therefore, the required number of states is the higher of the two requirements.")

smallest_coalition_size = max(states_for_proposal, states_for_ratification)

# Output the final calculation and result
print(f"Smallest number of states in a strong coalition = max({states_for_proposal}, {states_for_ratification})")
print(f"\nThe smallest number of states that could form a strong coalition is {smallest_coalition_size}.")
