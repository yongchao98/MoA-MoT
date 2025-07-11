# Define the parameters 'a' and 'b' from the curve equations x=cos(at) and y=sin(bt)
a = 9
b = 5

# For a Lissajous curve of the form x=cos(at), y=sin(bt) where a and b are
# coprime odd integers, the number of self-intersection points is given by
# the formula (a-1)*(b-1).
num_intersections = (a - 1) * (b - 1)

# Print the calculation and the final result.
# The request is to output each number in the final equation.
print(f"The curve is defined by (cos({a}t), sin({b}t)).")
print("This is a Lissajous curve with parameters a and b being coprime odd integers.")
print("The formula for the number of self-intersections in this case is (a-1) * (b-1).")
print(f"Calculation: ({a} - 1) * ({b} - 1) = {a-1} * {b-1} = {num_intersections}")
print(f"The number of self-intersection points is: {num_intersections}")
