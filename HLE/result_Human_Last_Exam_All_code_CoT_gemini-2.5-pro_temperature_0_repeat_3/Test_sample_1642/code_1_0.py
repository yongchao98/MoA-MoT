# Define the size of the chessboard.
board_size = 16

# Explain the logic based on the N-Queens problem.
print("The problem is to find the maximum number m for m white and m black queens")
print(f"on a {board_size}x{board_size} board where no queen attacks another.")
print("\nStep 1: Determine the total number of queens.")
print("The total number of queens is the sum of white and black queens: m + m = 2 * m.")

print("\nStep 2: Relate to the general N-Queens problem.")
print("A well-known result states that the maximum number of non-attacking queens")
print(f"you can place on a {board_size}x{board_size} board is {board_size}.")
max_total_queens = board_size

print("\nStep 3: Formulate the equation.")
print("We set the total number of queens (2 * m) equal to the maximum possible number of queens.")
# The equation is 2 * m = 16
m_coefficient = 2
print(f"The equation is: {m_coefficient} * m = {max_total_queens}")

print("\nStep 4: Solve for m.")
# Solve m = 16 / 2
final_m = max_total_queens // m_coefficient
print(f"m = {max_total_queens} / {m_coefficient}")
print(f"m = {final_m}")

print(f"\nTherefore, the maximum number m is {final_m}.")