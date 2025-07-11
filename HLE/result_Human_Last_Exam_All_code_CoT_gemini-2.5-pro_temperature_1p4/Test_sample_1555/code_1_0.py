import math

# This script calculates the leading-order asymptotic mass ratio between the lightest
# and next-to-lightest stable particles in the CP(N-1) model spectrum.

# In the large-N analysis of the CP(N-1) model, the spectrum contains a tower of
# stable particles, indexed by an integer k=1, 2, ... The mass of the k-th
# particle, M_k, is related to the mass of the lightest particle, M_1.

# The exact mass ratio is given by the formula:
# M_k / M_1 = sin(k * pi / N) / sin(pi / N)

# We are asked for the ratio of the subsequent higher excitation (k=2) to the
# lightest one (k=1).
print("The mass ratio R between the k=2 and k=1 states is R = M_2 / M_1.")

print("\nStep 1: State the formula for the mass ratio for k=2.")
print("R = sin(2 * pi / N) / sin(pi / N)")

print("\nStep 2: Simplify the expression using the trigonometric double-angle identity sin(2x) = 2*sin(x)*cos(x).")
print("Applying this identity with x = pi/N, the formula becomes:")
print("R = [2 * sin(pi / N) * cos(pi / N)] / sin(pi / N)")
print("R = 2 * cos(pi / N)")

print("\nStep 3: Determine the leading-order asymptotic value by taking the limit as N approaches infinity.")
print("We need to compute: lim (N -> inf) [2 * cos(pi / N)]")

print("\nStep 4: As N approaches infinity, the argument of cosine, (pi / N), approaches 0.")
print("The limit becomes: 2 * cos(0)")

print("\nStep 5: Calculate the final numerical result.")
# Define the constant from the equation
coefficient = 2
# Evaluate cos(0)
cos_of_zero = math.cos(0)
# Calculate the final result
final_ratio = coefficient * cos_of_zero

print(f"The value of cos(0) is {cos_of_zero}.")
print(f"The final equation is: {coefficient} * {cos_of_zero}")
print(f"The result is: {final_ratio}")

print("\nThus, the leading-order asymptotic mass ratio is 2.")
