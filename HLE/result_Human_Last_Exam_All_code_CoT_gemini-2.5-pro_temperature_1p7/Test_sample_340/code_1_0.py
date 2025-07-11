# Based on the analysis of the three-check chess position, a forced win for White can be achieved in a sequence of 3 moves.
# The plan is to deliver three consecutive checks that Black cannot parry indefinitely.
#
# The winning sequence is as follows:
# 1. White plays Qxf7+. This is Check #1. Black is forced to reply Kxf7.
# 2. White plays Bxd7+. This is Check #2. Black must respond.
# 3. No matter how Black responds, White can deliver a final check on the third move. For instance, if Black plays Qxd7, White plays Bd8+, which is Check #3.
#
# Each of these checking moves for White is counted as one move in the total.

# The move that delivers the first check.
moves_for_check_1 = 1

# The second move delivers the second check.
moves_for_check_2 = 1

# The third move delivers the final check.
moves_for_check_3 = 1

# The minimal number of moves for White to win is the sum of these three moves.
minimal_moves_to_win = moves_for_check_1 + moves_for_check_2 + moves_for_check_3

# The problem asks for the final answer as an integer, which is the result of this sequence.
# To satisfy the instruction "output each number in the final equation", the following line demonstrates the calculation.
print(f"The calculation for the minimal moves is: {moves_for_check_1} + {moves_for_check_2} + {moves_for_check_3} = {minimal_moves_to_win}")

# The final answer is the integer result.
print("\nThe minimal amount of moves by white to win is:")
print(minimal_moves_to_win)