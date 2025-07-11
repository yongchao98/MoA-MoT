# Define the parameters 'a' and 'b' from the curve equations x=cos(at) and y=sin(bt)
a = 9
b = 5

# The number of self-intersections for a Lissajous curve with coprime odd parameters a and b
# is calculated using the formula: (a - 1) * (b - 1) / 2.
# We first check if the parameters meet the criteria.
#
# 1. Check if a and b are integers (they are).
# 2. Check if a and b are odd.
is_a_odd = a % 2 != 0
is_b_odd = b % 2 != 0
#
# 3. Check if a and b are coprime. We can use a simple function for gcd.
def gcd(p, q):
    while q:
        p, q = q, p % q
    return p

are_coprime = gcd(a, b) == 1

# If all conditions are met, we can apply the formula.
if isinstance(a, int) and isinstance(b, int) and is_a_odd and is_b_odd and are_coprime:
    # Calculate the number of self-intersections
    intersections = (a - 1) * (b - 1) / 2

    # Print the equation with the numbers plugged in, and the final result.
    # The result should be an integer.
    print(f"The number of self-intersection points for the curve (cos({a}t), sin({b}t)) is calculated as:")
    print(f"({a} - 1) * ({b} - 1) / 2 = {int(intersections)}")
else:
    print("The formula for coprime odd integers does not apply to the given parameters.")
