import cmath

# The problem states the probability of selecting a reachable cell is 66%.
# P = 0.66
# Total cells = n*n
# Let R be the number of reachable cells. P = R / (n*n)

# Based on our model:
# Reachable white cells = (n*n)/2
# Reachable black cells = 2*n - 4
# So, R = (n*n)/2 + 2*n - 4

# The equation to solve is:
# ((n*n)/2 + 2*n - 4) / (n*n) = 0.66
# (n*n)/2 + 2*n - 4 = 0.66 * n*n
# 0.5*n*n + 2*n - 4 = 0.66 * n*n
# 0 = 0.16*n*n - 2*n + 4

# This is a quadratic equation of the form a*n^2 + b*n + c = 0
a = 0.16
b = -2
c = 4

# Calculate the discriminant
discriminant = (b**2) - 4*(a*c)

# Find the two solutions for n
sol1 = (-b - cmath.sqrt(discriminant)) / (2*a)
sol2 = (-b + cmath.sqrt(discriminant)) / (2*a)

# We are looking for an even integer solution.
# The real part of the solutions will give us the potential values for n.
n1 = sol1.real
n2 = sol2.real

final_n = 0
# Check which solution is a positive, even integer
if n1 > 0 and n1 == int(n1) and int(n1) % 2 == 0:
    final_n = int(n1)
elif n2 > 0 and n2 == int(n2) and int(n2) % 2 == 0:
    final_n = int(n2)

# Print the final equation with the determined value of n
# (n^2 / 2 + 2n - 4) / n^2 = 0.66
# We present it as 2n^2 - 25n + 50 = 0 for clarity of integer coefficients
# derived from 0.16n^2 - 2n + 4 = 0 by multiplying by 12.5 (or 100/8)
# Let's derive it cleanly:
# 0.16*n^2 - 2*n + 4 = 0
# multiply by 100: 16*n^2 - 200*n + 400 = 0
# divide by 8: 2*n^2 - 25*n + 50 = 0
final_a = 2
final_b = -25
final_c = 50

print(f"The model for the number of reachable cells is (n^2)/2 + 2n - 4.")
print(f"This leads to the quadratic equation for n: {final_a}n^2 + ({final_b})n + {final_c} = 0.")
print(f"Solving this equation gives two possible values for n: {n1:.2f} and {n2:.2f}.")
print(f"Since n must be an even integer, the only valid solution is {final_n}.")
print(f"Let's check the probability for n = {final_n}:")
reachable_cells = (final_n**2 / 2) + 2*final_n - 4
total_cells = final_n**2
probability = reachable_cells / total_cells
print(f"Reachable cells = ({final_n}^2 / 2) + 2*{final_n} - 4 = {int(reachable_cells)}")
print(f"Total cells = {final_n}^2 = {total_cells}")
print(f"Probability = {int(reachable_cells)} / {total_cells} = {probability}")

<<<10>>>