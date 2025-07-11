# Define the number of turns and their respective rotations.
fouette_turns_count = 8
fouette_rotations_per_turn = 1

pirouette_count = 1 
pirouette_rotations_per_turn = 3 # This is a triple pirouette

chaine_turns_count = 6
chaine_rotations_per_turn = 1

# Calculate the total rotations for each type of turn.
total_fouette_rotations = fouette_turns_count * fouette_rotations_per_turn
total_pirouette_rotations = pirouette_count * pirouette_rotations_per_turn
total_chaine_rotations = chaine_turns_count * chaine_rotations_per_turn

# Calculate the grand total number of rotations.
grand_total_rotations = total_fouette_rotations + total_pirouette_rotations + total_chaine_rotations

# Print the final equation as requested, showing each component.
print(f"The total number of rotations is the sum of rotations from fouettés, pirouettes, and chaînés.")
print(f"{total_fouette_rotations} + {total_pirouette_rotations} + {total_chaine_rotations} = {grand_total_rotations}")