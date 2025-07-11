# In triangle ABC, side lengths are denoted by a, b, and c.
# a is the length of side BC, opposite to angle A.
# b is the length of side CA, opposite to angle B.
# c is the length of side AB, opposite to angle C.

# Based on the geometric derivation, the ratio BM/MI is (a + c) / b.

# The final equation for the ratio is:
numerator_term_1 = 'a'
numerator_term_2 = 'c'
operator = '+'
denominator_term = 'b'

print("The formula for the ratio BM/MI is expressed in terms of the side lengths a, b, and c.")
print("Final Equation:")
print(f"Ratio = ({numerator_term_1} {operator} {numerator_term_2}) / {denominator_term}")

# This can also be defined as a function for calculation with specific values
def get_ratio(a, b, c):
  """
  Calculates the ratio BM/MI for a triangle with side lengths a, b, c.
  Checks for the triangle inequality.
  """
  if a + b > c and a + c > b and b + c > a:
    return (a + c) / b
  else:
    return "The given side lengths do not form a valid triangle."

# Example usage with a triangle of sides a=3, b=4, c=5 (a right triangle)
# ratio_value = get_ratio(3, 4, 5)
# print(f"\nFor an example triangle with a=3, b=4, c=5:")
# print(f"The ratio BM/MI is ({3} + {5}) / {4} = {ratio_value}")
