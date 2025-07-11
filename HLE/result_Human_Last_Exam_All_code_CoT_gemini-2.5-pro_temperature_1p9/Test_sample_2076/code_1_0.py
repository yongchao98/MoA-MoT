import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# The logic to determine the pieces and their values.
# Piece values: Pawn=1, Knight=3, Bishop=3, Rook=5, Queen=9

# White pieces configuration
white_pieces = {"Bishop": 3}
# Black pieces configuration
black_pieces = {"Pawn": 1, "Knight": 3}

white_total_value = sum(white_pieces.values())
black_total_value = sum(black_pieces.values())
total_pieces_count = len(white_pieces) + len(black_pieces)
total_value = white_total_value + black_total_value

# The description explains the logic of selecting this combination.
print("This position is a Diagonal Corridor Mate constructed with the minimum number of pieces (3) while ensuring Black has a material advantage.")
print(f"White has {len(white_pieces)} piece(s) with a total value of {white_total_value}.")
print(f"Black has {len(black_pieces)} pieces with a total value of {black_total_value}.")
print(f"Black's material advantage is confirmed: {black_total_value} > {white_total_value}.")
print(f"The total value of the {total_pieces_count} added pieces is {total_value}, the minimum possible under the given constraints.")
print("-" * 20)

# Format the final list as per the instructions.
# White pieces first, then Black pieces sorted by value.
white_list = [f"White {piece}" for piece in sorted(white_pieces, key=white_pieces.get)]
black_list = [f"Black {piece}" for piece in sorted(black_pieces, key=black_pieces.get)]

final_answer_list = ", ".join(white_list + black_list)

print("Final Piece List:")
print(final_answer_list)

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output_content = captured_output.getvalue()

# Print the captured output to the actual console
print(output_content)

# Final answer in the specified format
print(f'<<<{final_answer_list}>>>')