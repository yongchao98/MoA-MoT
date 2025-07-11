# The best sequence is a smothered mate, which is the fastest forced checkmate in this position.
# We will construct the string representing the move sequence and print it.

# Define the parts of the algebraic notation sequence
black_turn_1 = "1..."
black_move_1 = "Qg1+"
white_turn_2 = "2."
white_move_2 = "Rxg1"
black_move_2 = "Nf2#"

# Combine the parts into the final equation string.
# This ensures we output each number and symbol in the final equation as requested.
final_equation = f"{black_turn_1} {black_move_1} {white_turn_2} {white_move_2} {black_move_2}"

# Print the final sequence
print(final_equation)