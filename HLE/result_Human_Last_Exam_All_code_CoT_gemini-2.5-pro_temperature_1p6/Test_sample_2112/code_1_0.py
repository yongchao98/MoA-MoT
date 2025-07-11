import math

# The problem reduces to solving for r_0 in one of several cases derived from the algebraic structure.
# The case that yields a solution r_0 > 15 is when T(r_0) = 1/sqrt(2), where T(r) = (3r-37)/(r+4).
# So we need to solve the equation: (3*r_0 - 37) / (r_0 + 4) = 1 / sqrt(2)
# Rearranging this gives: sqrt(2) * (3*r_0 - 37) = r_0 + 4
# r_0 * (3*sqrt(2) - 1) = 4 + 37*sqrt(2)
# r_0 = (4 + 37*sqrt(2)) / (3*sqrt(2) - 1)

# Let's calculate the terms of the rearranged linear equation for r_0: A*r_0 = B
A = 3 * math.sqrt(2) - 1
B = 4 + 37 * math.sqrt(2)

# Solve for r_0
r_0 = B / A

# The problem asks to output the numbers in the final equation.
# We will present the linear equation from which r_0 is solved.
# The equation is (3*sqrt(2)-1)*r_0 - (4+37*sqrt(2)) = 0
term1_coeff = A
term2_const = -B

# To fulfill the requirement of showing each number, we will express r_0 as a rationalized expression.
# r_0 = (4 + 37*sqrt(2))*(3*sqrt(2)+1) / ((3*sqrt(2)-1)*(3*sqrt(2)+1))
# Denominator = (3*sqrt(2))^2 - 1^2 = 9*2 - 1 = 17
# Numerator = 4*3*sqrt(2) + 4*1 + 37*sqrt(2)*3*sqrt(2) + 37*sqrt(2)*1
# Numerator = 12*sqrt(2) + 4 + 37*3*2 + 37*sqrt(2) = 12*sqrt(2) + 4 + 222 + 37*sqrt(2) = 226 + 49*sqrt(2)
# So r_0 = (226 + 49*sqrt(2)) / 17

numerator_const = 226
numerator_sqrt_coeff = 49
denominator = 17

print("The problem simplifies to solving an equation for r_0.")
print("The specific case that yields a valid solution gives the equation T(r_0) = 1/sqrt(2), where T(r) = (3r-37)/(r+4).")
print("Rearranging this leads to a linear equation for r_0 of the form A*r_0 = B.")
print(f"The equation is: ({3}*sqrt({2}) - {1}) * r_0 = {4} + {37}*sqrt({2})")
print(f"This can be written as {A} * r_0 = {B}")
print("Solving for r_0 and rationalizing the denominator gives r_0 = (c + d*sqrt(2))/e")
print(f"The equation with integer coefficients for the rationalized form is:")
print(f"{denominator}*r_0 - ({numerator_const} + {numerator_sqrt_coeff}*sqrt({2})) = 0")
print("\nThe value of the radial distance r_0 is:")
print(r_0)
