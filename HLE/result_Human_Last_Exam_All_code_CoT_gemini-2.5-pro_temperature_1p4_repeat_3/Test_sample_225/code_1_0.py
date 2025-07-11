# Define the parameters 'a' and 'b' from the curve equations x=cos(at) and y=sin(bt)
a = 9
b = 5

# The number of self-intersection points for a Lissajous curve with coprime integer
# parameters a and b is given by the formula (a-1)(b-1)/2.
# We first check if a and b are coprime. The greatest common divisor (GCD) can be found using the math module.
import math
if math.gcd(a, b) == 1:
    # Calculate the intermediate values
    val_a = a - 1
    val_b = b - 1
    numerator = val_a * val_b
    
    # Calculate the final number of self-intersection points
    num_intersections = numerator / 2
    
    # Print the equation and the final answer as requested
    print(f"The number of self-intersection points is calculated by the formula (a-1)*(b-1)/2.")
    print(f"For a={a} and b={b}:")
    print(f"({a} - 1) * ({b} - 1) / 2 = {int(num_intersections)}")
else:
    print("The formula does not apply as 'a' and 'b' are not coprime.")
