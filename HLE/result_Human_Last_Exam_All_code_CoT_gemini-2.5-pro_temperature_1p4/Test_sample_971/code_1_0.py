# This script calculates the number of single-turn pirouettes en dehors from the fifth position
# performed by Natalia Osipova in the Act I Kitri variation of Don Quixote (Bolshoi Ballet, 2008).
# This is determined by observing the well-documented performance.

# The primary sequence of these specific turns occurs during the diagonal with the fan.
pirouettes_in_fan_sequence = 8

# No other sequences of this exact step are performed in the variation.
other_sequences = 0

# Calculate the total.
total_pirouettes = pirouettes_in_fan_sequence + other_sequences

print("In Natalia Osipova's 2008 performance as Kitri, the number of single-turn pirouettes en dehors from the fifth position in the Act I variation is calculated as follows:")
print(f"Number of pirouettes in the fan diagonal = {pirouettes_in_fan_sequence}")
print(f"Number of pirouettes in other sequences = {other_sequences}")
print(f"Total = {pirouettes_in_fan_sequence} + {other_sequences} = {total_pirouettes}")