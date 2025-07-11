# The solution to the system of differential equations blows up under certain
# initial conditions. Given x(0) > 1, the condition on y(0) is determined
# by a special curve in the phase plane known as a separatrix.

# Through analytical methods, the equation for this separatrix has been
# found to be: y**2 = 2*x + 1 - 3*x**(2/3).

# A solution blows up if the initial state (x(0), y(0)) is below this curve.
# This gives a clear condition on the initial value y(0).

# The numbers defining the separatrix equation are extracted here:
coeff_x = 2
constant_term = 1
coeff_power_term = 3
exponent_numerator = 2
exponent_denominator = 3

# The following code prints the blow-up condition using these numbers.
print("Given x(0) > 1, the solution of the system blows up if y(0) satisfies the following inequality:")

# We construct and print the final equation for the condition on y(0).
# The expression represents the upper boundary for y(0) for a given x(0).
condition_str = (
    f"y(0) < sqrt({coeff_x}*x(0) + {constant_term} - "
    f"{coeff_power_term}*x(0)**({exponent_numerator}/{exponent_denominator}))"
)

print(condition_str)