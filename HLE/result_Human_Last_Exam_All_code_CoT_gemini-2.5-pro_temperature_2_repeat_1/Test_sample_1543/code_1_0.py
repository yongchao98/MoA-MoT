# The FEN describes a Capablanca Chess position:
# 9k/5c1pb1/10/10/10/3Q6/PP5A2/K9 w - - 0 1
# White has a forced checkmate sequence.

# Move 1 by White: Archbishop from h2 to j3 (as a knight)
white_move_1 = "1. Aj3"

# After 1. Aj3, Black is faced with unstoppable threats.
# Let's consider Black's best moves and White's mating response.

# Scenario 1: Black moves the King
# If Black plays 1... Kh8, White delivers mate with 2. Qf5#
# If Black plays 1... Ki8, White delivers mate with 2. Qxh7#
# In this scenario, White takes 2 moves to win.
winning_moves_scenario_1 = 2

# Scenario 2: Black defends with the Chancellor
# If Black plays 1... Cg7 (to prevent a Qg7 mate), White plays 2. Qxg7#
# In this scenario, White also takes 2 moves to win.
winning_moves_scenario_2 = 2

# In all optimal lines, White delivers checkmate on the second move.
# Therefore, the minimal number of moves for White to win is 2.

# Representing this as an equation:
# (1st move by White) + (2nd move by White) = 2 total moves by White
move1 = 1
move2 = 1
total_moves = move1 + move2

print(f"White's first move: {white_move_1}")
print(f"Black must respond. For any optimal response from Black, White delivers checkmate on the next move.")
print(f"For example, if Black plays 1... Ki8, White wins with 2. Qxh7#.")
print("The minimal number of moves by White to win is therefore 2.")
print(f"Final equation: {move1} + {move2} = {total_moves}")

# The final answer is the integer representing the number of moves.
print(f"So the answer is: {total_moves}")
