# Natalia Osipova as Kitri in Don Quixote, Act I Variation (Bolshoi Ballet, 2008)
# This script calculates the number of single-turn pirouettes en dehors
# from the fifth position performed consecutively in her famous sequence.

# Each turn is counted as 1.
turn_1 = 1
turn_2 = 1
turn_3 = 1
turn_4 = 1
turn_5 = 1
turn_6 = 1
turn_7 = 1
turn_8 = 1

# We create a list to hold each turn for the final printed equation.
all_turns = [turn_1, turn_2, turn_3, turn_4, turn_5, turn_6, turn_7, turn_8]

# Calculate the total number of turns.
total_turns = sum(all_turns)

# Create the equation string.
equation = " + ".join(str(turn) for turn in all_turns)

# Print the final result in the format of an equation.
print(f"The number of single-turn pirouettes performed is:")
print(f"{equation} = {total_turns}")