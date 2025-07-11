# Side lengths of the triangle ABC
# a is the length of side BC
# b is the length of side AC
# c is the length of side AB

# We can use example values for the side lengths to demonstrate.
# Let's use a common 3-4-5 right triangle for this example.
a = 3.0
b = 4.0
c = 5.0

# The problem asks for the ratio BM/MI.
# The step-by-step derivation shows that this ratio is equal to (a + c) / b.

# Check if the side lengths form a valid triangle
if a + b > c and a + c > b and b + c > a:
    # Symbolic representation of the formula
    numerator_symbolic = "a + c"
    denominator_symbolic = "b"
    
    # Calculate the numerical value of the ratio
    ratio_value = (a + c) / b
    
    # Print the final derived expression
    print(f"The ratio BM/MI can be expressed in terms of the side lengths a, b, and c as:")
    print(f"({numerator_symbolic}) / ({denominator_symbolic})")
    
    # Print the final equation with the example values
    print(f"\nFor a triangle with side lengths a = {a}, b = {b}, c = {c}:")
    print("The final equation is:")
    print(f"BM / MI = ({a} + {c}) / {b} = {ratio_value}")
else:
    print("The given side lengths do not form a valid triangle.")
