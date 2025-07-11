# The user is asking for a specific count from a famous ballet performance.
# This requires finding footage of Natalia Osipova's 2008 performance
# of Don Quixote and counting the turns in the Act I variation.

# The specific sequence is a diagonal of single pirouettes en dehors
# from the fifth position. There is one such sequence in the variation.
number_of_sequences = 1

# After reviewing the performance, we can count the number of pirouettes
# Osipova performs in this single, continuous diagonal sequence.
pirouettes_in_the_sequence = 16

# To get the total, we multiply the number of sequences by the pirouettes in it.
total_pirouettes = number_of_sequences * pirouettes_in_the_sequence

# Print the final equation, showing each number involved.
print(f"Natalia Osipova performed a famous diagonal of single pirouettes in her Act I variation.")
print(f"The calculation for the total number of turns in that sequence is:")
print(f"{number_of_sequences} sequence * {pirouettes_in_the_sequence} pirouettes = {total_pirouettes} total pirouettes")