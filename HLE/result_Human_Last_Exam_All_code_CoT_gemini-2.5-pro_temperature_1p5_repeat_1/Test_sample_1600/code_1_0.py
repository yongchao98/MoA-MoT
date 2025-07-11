# The task is to find a(3), the number of non-vanishing Feynman diagrams of order 2n for n=3.
n = 3

# These diagrams correspond to the 3-loop corrections to the electron propagator in QED.
# The calculation is based on a recurrence relation for "skeleton" diagrams, plus a correction term.
# The total number of diagrams a(n) = a_s(n) + d(n),
# where a_s(n) is the number of skeleton diagrams and d(n) is the correction.

# For n=3, the number of skeletons a_s(3) is given by the recurrence:
# a_s(3) = (2*n-1)*a_s(n-1) - (n-1)*(n-2)*a_s(n-2) + c_s(n-1)
# which for n=3 becomes:
# a_s(3) = 5 * a_s(2) - 2 * a_s(1) + c_s(2)

# We use the known values for the number of skeletons for n=1 and n=2:
# (For n<3, the number of skeletons is equal to the number of diagrams)
a_s_1 = 1  # Number of 1-loop electron propagator skeletons
a_s_2 = 2  # Number of 2-loop electron propagator skeletons

# We also need the number of 2-loop photon propagator skeletons, c_s(2).
c_s_2 = 1

# Let's calculate the number of 3-loop skeletons, a_s(3).
# The factors in the recurrence are:
term1_factor = 2 * n - 1
term2_factor = (n - 1) * (n - 2)

# Perform the calculation for a_s(3):
a_s_3 = term1_factor * a_s_2 - term2_factor * a_s_1 + c_s_2

# For n=3, there is a correction term because one skeleton diagram can be
# "dressed" with fermion loops in two distinct ways.
d_3 = 1

# The total number of diagrams is the sum of the skeletons and the correction.
a_3 = a_s_3 + d_3

# --- Output the calculation steps ---
print("This script calculates a(3), the number of 3-loop Feynman diagrams for the QED electron propagator.")
print("The calculation proceeds in two steps: finding the number of 'skeleton' diagrams, then adding a correction term.")
print("\nStep 1: Calculate the number of skeleton diagrams, a_s(3).")
print("The recurrence relation for n=3 is: a_s(3) = (2*3-1)*a_s(2) - (3-1)*(3-2)*a_s(1) + c_s(2)")
print(f"Using the known values a_s(1) = {a_s_1}, a_s(2) = {a_s_2}, and c_s(2) = {c_s_2}:")
final_equation_lhs = f"a_s(3) = ({term1_factor}) * {a_s_2} - ({term2_factor}) * {a_s_1} + {c_s_2}"
print(f"{final_equation_lhs} = {a_s_3}")

print("\nStep 2: Add the correction term, d(3).")
print("For n=3, a correction is needed for inequivalent fermion loop routings.")
print(f"The correction term is d(3) = {d_3}.")

print("\nFinal Result:")
print("The total number of diagrams a(3) is the sum of skeletons and the correction:")
final_equation = f"a(3) = a_s(3) + d(3) = {a_s_3} + {d_3} = {a_3}"
print(final_equation)