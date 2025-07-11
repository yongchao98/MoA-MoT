# This problem is a chess puzzle that requires logical deduction rather than a simple calculation.
# The goal is to find the minimum number of moves for White to win a Three-Check chess game
# from the given position, assuming optimal play from Black.

# The winning sequence of moves for White is as follows:
# 1. Qxf7+ (Check #1). Black must play Kxf7.
# 2. Bc4+ (Check #2). Black's optimal defense to prolong the game is Kf8.
# 3. O-O-O (A setup move, as there's no immediate 3rd check). Black plays a defensive move.
# 4. Rxd7+ (Check #3). White wins.

# Therefore, the minimal number of moves by White to win is 4.
minimal_white_moves_to_win = 4

# The problem asks for a single non-negative integer as the answer.
print(minimal_white_moves_to_win)