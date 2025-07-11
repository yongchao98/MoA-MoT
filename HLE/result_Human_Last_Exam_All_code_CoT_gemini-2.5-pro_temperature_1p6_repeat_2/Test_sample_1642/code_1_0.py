# Define the size of the chessboard
N = 16

# The total number of queens is m white queens + m black queens = 2*m.
# The maximum number of non-attacking queens on an N x N board is N.
# To find the maximum value for m, we set the total number of queens equal to N.
# The equation is: 2 * m = N

# Solve for m
m = N / 2

# Ensure m is an integer for the final output
m_int = int(m)

# Print the final result including the equation with the solved value of m.
print(f"The maximum number of queens (m) of each color is found by solving the equation 2 * m = {N}.")
print(f"Solving for m gives m = {N} / 2 = {m_int}.")
print("The final equation with the numbers is:")
print(f"2 * {m_int} = {N}")
