# This problem is a case of the Furstenberg set problem from geometric measure theory.
# Let s be the required minimal dimension of the intersection of the set C with a line.
# In the problem statement, this value is 1/2.
s = 0.5

# The minimal possible dimension of the set C, let's call it d(s), is conjectured to follow
# a two-part formula depending on the value of s.
# This conjecture was formulated by Thomas Wolff.
# The lower bound d(s) >= max(2*s, s + 1/2) is proven, and for s=0.5, this bound is 1.
# The conjecture matches this lower bound at s=0.5.

# According to the conjecture:
# If s is between 0 and 1/2 (inclusive), the dimension is 2*s.
# If s is between 1/2 and 1 (inclusive), the dimension is s + 1/2.
if 0 <= s <= 0.5:
  minimal_dimension = 2 * s
elif 0.5 < s <= 1:
  minimal_dimension = s + 0.5
else:
  minimal_dimension = -1 # s is out of the standard range for this problem

# We print the parameters of the equation used and the result.
# The formula used here is d(s) = 2 * s, since s = 0.5.
factor = 2
print(factor)
print(s)
print(minimal_dimension)
