# The user wants to find the correct move for Black on a Go board.
# The goal is to make a move that ensures a Black group survives AND captures two White stones.
#
# Board analysis:
# 1. The endangered Black group is at coordinates (0, 2) and (1, 2). It has only two liberties ((0, 3) and (1, 3)) and can be captured.
# 2. The target White stones are at (0, 1) and (1, 1).
# 3. The solution is a classic Go tactic called a "throw-in" sacrifice, leading to a "snapback".
#
# The move sequence:
# 1. Black plays at (0, 0). This is a sacrifice.
# 2. White must capture the Black stone at (0, 0) by playing at (1, 0). If White doesn't, Black will play at (1, 0) and put the White group in atari (a state of having only one liberty).
# 3. After White plays at (1, 0), the White group at (0, 1), (1, 1) and the newly captured area are now in atari themselves. Their single, shared liberty is at (2, 1).
# 4. It's Black's turn again. Black plays at (2, 1) and captures the entire White group.
# 5. By capturing the White stones, the originally endangered Black group is now safe and connected.
#
# Therefore, the correct move for Black to start this sequence is (0, 0).
# The coordinate system is (row, column) starting from the top-left at (0, 0).

print("(0, 0)")