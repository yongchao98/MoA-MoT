# Step 1: Define the given number of undisturbed curves from the poem's scenario.
undisturbed_curves = 90

# Step 2: Set up the equation to find the original total number of curves (T).
# T = undisturbed + (3/8)*T + (1/4)*T
# T = 90 + (3/8 + 2/8)*T
# T = 90 + (5/8)*T
# T - (5/8)*T = 90
# (3/8)*T = 90
original_total = (90 * 8) // 3

# Step 3: Calculate the number of curves that "lost shape". This is the group
# that needs to be redrawn, as it's the only one that works with the
# subsequent fractions (1/5 and 2/9), as its size must be a multiple of 45.
drawn_curves = (3 * original_total) // 8

# Let's verify this number is a multiple of 5 and 9.
# Parabolic curves = 1/5 * drawn_curves
# Elliptical curves = 2/9 * drawn_curves
# Both must be integers.
parabolic_curves = drawn_curves // 5
elliptical_curves = (2 * drawn_curves) // 9
circular_curves = drawn_curves - parabolic_curves - elliptical_curves

# Step 4: Calculate the final total number of curves in the restored pattern.
# This is the sum of the undisturbed curves and the newly drawn curves.
final_total = undisturbed_curves + drawn_curves

# Step 5: Print the final answer, showing the breakdown.
# The question asks for the total curves. The code will output the final equation.
print(f"The number of undisturbed curves is {undisturbed_curves}.")
print(f"The number of curves the master must draw is {drawn_curves}.")
print("The final pattern's total number of curves is the sum of the undisturbed curves and the newly drawn curves.")
print(f"{undisturbed_curves} + {drawn_curves} = {final_total}")