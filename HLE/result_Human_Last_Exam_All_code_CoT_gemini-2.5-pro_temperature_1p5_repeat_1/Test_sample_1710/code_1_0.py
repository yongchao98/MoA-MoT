# The user wants to find the correct move for Black on a Go board.
# The objective is for Black to survive and capture two specific White stones.
#
# 1. Identify the target: The two isolated White stones are at coordinates (2, 6) and (3, 6).
# 2. Analyze the situation: Capturing these stones will create a large, secure "eye" for the
#    surrounding Black group, ensuring it survives. The capture requires a "geta" (net) technique,
#    as a direct atari is not possible.
# 3. Identify candidate moves: The key escape points for the White stones are (1, 6), (2, 5), (3, 5), and (4, 6).
#    The candidate moves to form a net are typically knight's moves away, which are (1, 5) and (4, 5).
# 4. Evaluate candidates:
#    - A move at (1, 5) successfully cuts off all of White's escape paths. If White tries to run up or left,
#      Black can block. If White tries to run down, Black can play a follow-up move at (4, 5) to complete the net.
#    - A move at (4, 5) also successfully cuts off all escape paths. If White tries to run down, Black can
#      block. If White tries to run up, Black has follow-up moves (like (2,4)) to cut White's connection to other
#      groups and secure the capture.
# 5. Conclusion: Both (1, 5) and (4, 5) are correct moves.
# 6. Formatting: The answers must be ordered lexicographically. (1, 5) comes before (4, 5).

# The final answer consists of two coordinate pairs.
answer_1_row = 1
answer_1_col = 5
answer_2_row = 4
answer_2_col = 5

# Print the final answer in the specified format: (row, col), (row, col)
print(f"({answer_1_row}, {answer_1_col}), ({answer_2_row}, {answer_2_col})")