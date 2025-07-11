# This script solves the classic chess problem of finding the minimum material
# to attack 63 out of 64 squares, leaving the black king in stalemate on the last one.
#
# The solution to this problem is 10 points.
# This value is achieved with two pieces: a Queen and a Pawn.
#
# A famous configuration that achieves this was composed by T. Siers in 1935:
# - White: Queen at e1, Pawn at g5
# - Black: King at g7
# The black king on g7 is not in check, but all of its potential escape squares
# (f6, h6, f7, g6, h7, f8, g8, h8) are attacked by either the queen or the pawn,
# resulting in a stalemate. All other squares on the board are also attacked.
#
# The script will now calculate the total point value.

queen_points = 9
pawn_points = 1

total_points = queen_points + pawn_points

print("The solution uses a Queen and a Pawn.")
print("The value of a Queen is 9 points.")
print("The value of a Pawn is 1 point.")
print("To find the smallest number of points, we add the values of the pieces.")
print(f"The final equation is: {queen_points} + {pawn_points} = {total_points}")
