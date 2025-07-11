import math

# The parameters for the Lissajous curve x = cos(at), y = sin(bt)
a = 9
b = 5

# The number of self-intersection points for a Lissajous curve where a and b are
# coprime integers is given by the formula (a - 1)(b - 1) / 2.
# We first confirm that a and b are coprime.
if math.gcd(a, b) == 1:
    # Calculate the individual terms of the formula
    term1 = a - 1
    term2 = b - 1
    numerator = term1 * term2
    num_intersections = numerator // 2

    # Print the step-by-step calculation as a final equation
    print("The number of self-intersection points is calculated using the formula (a-1)(b-1)/2:")
    print(f"({a} - 1) * ({b} - 1) / 2 = {term1} * {term2} / 2 = {numerator} / 2 = {num_intersections}")
else:
    # This part would execute if a and b were not coprime.
    print(f"The formula for coprime parameters does not apply because gcd({a}, {b}) is not 1.")
