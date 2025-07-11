# The problem asks for the range of y(0) given x(0) > 1 that causes the solution to blow up.
# Based on the analytical derivation, the condition for blow-up is an inequality involving the initial conditions.
# This script prints the resulting inequality.

# The numbers in the final inequality are:
zero = 0
one = 1
two = 2
three = 3

# The condition is y(0) < sqrt(2*x(0) + 1 - 3*x(0)^(2/3)).
# We print this final inequality for y(0), showing each number explicitly as requested.
print(f"Given x({zero}) > {one}, the solution to the system of differential equations blows up if and only if y({zero}) satisfies the following inequality:")
print(f"y({zero}) < ({two}*x({zero}) + {one} - {three}*x({zero})**({two}/{three}))**({one}/{two})")