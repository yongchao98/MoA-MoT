# Define the parameters of the Lissajous curve
a = 9
b = 5

# The formula for the number of self-intersections for a curve
# (cos(at), sin(bt)) with coprime a and b is (a-1)(b-1)/2.

# Calculate the intermediate values for the formula
a_minus_1 = a - 1
b_minus_1 = b - 1
numerator = a_minus_1 * b_minus_1
num_intersections = numerator // 2

# Print the explanation and the step-by-step calculation
print(f"The curve is defined by the parametric equations (cos({a}t), sin({b}t)).")
print("To find the number of self-intersections, we use the formula N = (a - 1) * (b - 1) / 2.")
print(f"N = ({a} - 1) * ({b} - 1) / 2")
print(f"N = {a_minus_1} * {b_minus_1} / 2")
print(f"N = {numerator} / 2")
print(f"N = {num_intersections}")