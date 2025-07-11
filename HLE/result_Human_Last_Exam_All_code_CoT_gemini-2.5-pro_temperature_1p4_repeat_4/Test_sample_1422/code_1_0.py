# Set the parameters based on the problem statement
# L is the number of loops
L = 2
# n is the number of lines per vertex. We consider the common phi^4 theory.
n = 4
# E is the number of external lines. To minimize vertices, we choose the minimum E.
E = 0

print(f"To find the minimum number of vertices (V) in a two-loop Feynman diagram, we use the formula:")
print(f"V = (2*L - 2 + E) / (n - 2)")
print(f"We are given L = {L} (two-loop).")
print(f"For the simplest interacting scalar theory (phi^n), we choose n = {n}.")
print(f"To minimize V, we choose the minimum number of external lines, E = {E}.")
print("")

# Calculate the number of vertices V
numerator = 2 * L - 2 + E
denominator = n - 2
V = numerator / denominator

print("Plugging in the values, the equation is:")
print(f"V = (2 * {L} - 2 + {E}) / ({n} - 2)")
print(f"V = {numerator} / {denominator}")

# The final result must be an integer
min_vertices = int(V)

print(f"The minimum number of vertices is: {min_vertices}")
print("\nThis corresponds to the 'figure-eight' vacuum diagram in phi^4 theory.")

<<<1>>>