# The curve is described by the parametric equations x(t) = cos(9t) and y(t) = sin(5t).
# This is a type of curve called a Lissajous curve.
# For a Lissajous curve of the form x(t) = cos(at), y(t) = sin(bt), the number of
# self-intersection points is given by the formula (a-1)(b-1), provided that a and b
# are coprime integers (their greatest common divisor is 1).

# Define the parameters 'a' and 'b' from the equations.
a = 9
b = 5

# First, let's confirm that a and b are coprime. We can do this by finding their greatest common divisor (GCD).
import math
gcd = math.gcd(a, b)
if gcd != 1:
    print(f"The parameters a={a} and b={b} are not coprime (GCD={gcd}), so the formula may not apply.")
else:
    # Calculate the number of self-intersections using the formula.
    num_intersections = (a - 1) * (b - 1)

    # Print the explanation and the final equation with all its numbers.
    print(f"The parameters of the curve are a = {a} and b = {b}.")
    print(f"Since gcd({a}, {b}) = {gcd}, they are coprime and the formula (a-1)*(b-1) can be used.")
    print(f"The final calculation is: ({a} - 1) * ({b} - 1) = {a - 1} * {b - 1} = {num_intersections}")
    print(f"\nThe number of self-intersection points is {num_intersections}.")
