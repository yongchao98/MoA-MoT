import sympy

# Define symbols
x, t, c, i = sympy.symbols('x t c i')

# Step 1: Define the initial condition v(x,0) as f(x)
# f(x) = (1/(2*sqrt(6))) * (1 + exp(-x))**2
# We derived this from the problem statement.
sqrt6 = sympy.sqrt(6)
f = (1 / (2 * sqrt6)) * (1 + sympy.exp(-x))**2
print(f"The function f(x), which is v(x,0), is: {f}")

# Step 2: The result of the fractional derivative operation on f(x-ct)
# evaluated at x=ct is -i*sqrt(c)*f'(0).
# Here c = 6.
c_val = 6
print(f"The velocity of the travelling wave is c = {c_val}")

# Step 3: Differentiate f(x) with respect to x
dfdx = sympy.diff(f, x)
print(f"The derivative of f(x) is f'(x) = {dfdx}")

# Step 4: Evaluate the derivative at x=0
df_at_0 = dfdx.subs(x, 0)
print(f"The derivative at x=0 is f'(0) = {df_at_0}")
# Simplify the f'(0) expression to show the numbers in the equation
f_prime_0_val_numerator = -1 * (1 + sympy.exp(0)) * sympy.exp(0)
f_prime_0_val_denominator = sympy.sqrt(6)
print(f"So, f'(0) = ({f_prime_0_val_numerator}) / ({f_prime_0_val_denominator}) = {df_at_0}")


# Step 5: Calculate the final quantity Q = -i*sqrt(c)*f'(0)
i = sympy.I
Q = -i * sympy.sqrt(c_val) * df_at_0
print(f"\nThe final quantity Q is -i * sqrt({c_val}) * f'(0)")
# We can print out the final equation with the numbers plugged in.
print(f"Q = -i * ({sympy.sqrt(c_val)}) * ({df_at_0})")
final_answer = sympy.simplify(Q)
print(f"Q = {final_answer}")
