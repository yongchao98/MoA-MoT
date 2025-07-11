# The puzzle is to find the minimal number of moves for White to win.
# The analysis reveals a forced checkmate in 3 moves.
# The winning sequence is:
# 1. Qg6 hxg6
# 2. Af4+ Kh8
# 3. Ag6#

# The number of moves for White to win is 3.
# The final move is the 3rd move by white.

white_move_1 = "Qg6"
black_move_1 = "hxg6"
white_move_2 = "Af4+"
black_move_2 = "Kh8"
white_move_3 = "Ag6#"

number_of_white_moves_to_win = 3

print(f"The minimal amount of moves by White to win is {number_of_white_moves_to_win}.")
print(f"The winning line is: 1. {white_move_1} {black_move_1} 2. {white_move_2} {black_move_2} 3. {white_move_3}")
