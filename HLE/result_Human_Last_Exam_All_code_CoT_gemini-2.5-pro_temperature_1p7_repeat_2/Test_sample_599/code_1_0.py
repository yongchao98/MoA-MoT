# The problem is to find the 50th segmented number.
# A segmented number is a positive integer that cannot be expressed as the sum of
# two or more consecutive positive integers.
#
# Through mathematical analysis, it has been established that these numbers are
# precisely the powers of 2 (i.e., numbers of the form 2^p).
#
# The sequence of segmented numbers is:
# 1st element: 2^(1-1) = 2^0 = 1
# 2nd element: 2^(2-1) = 2^1 = 2
# 3rd element: 2^(3-1) = 2^2 = 4
# ...and so on.
#
# To find the 50th element, we use this pattern.
term_number = 50

# The formula for the Nth element is 2^(N-1).
power = term_number - 1
result = 2**power

# The final equation is result = 2^49.
# We print all the numbers involved in this equation.
print(f"The 50th segmented number is the result of the expression: 2^{power}")
print(f"The calculation is: 2^{power} = {result}")
