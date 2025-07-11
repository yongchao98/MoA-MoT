# Define the parameters 'a' and 'b' from the parametric equations
# x(t) = cos(a*t), y(t) = sin(b*t)
a = 9
b = 5

# For a Lissajous curve where 'a' and 'b' are coprime odd integers,
# the number of self-intersections is given by the formula (a-1)*(b-1)/2.

# Calculate the intermediate values for the final equation
val_a_minus_1 = a - 1
val_b_minus_1 = b - 1
numerator = val_a_minus_1 * val_b_minus_1
num_intersections = int(numerator / 2)

# Print the final equation with all the numbers, as requested
print("The number of self-intersection points is calculated using the formula N = (a-1)*(b-1)/2.")
print("Calculation:")
print(f"({a} - 1) * ({b} - 1) / 2 = {val_a_minus_1} * {val_b_minus_1} / 2 = {numerator} / 2 = {num_intersections}")

# Print the final answer clearly
print(f"\nThe number of self-intersection points is: {num_intersections}")