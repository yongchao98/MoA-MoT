import numpy as np

# Step 1: Find the coefficients of f(x)
print("Step 1: Determine the function f(x) = a*e^(2x) + b*e^x + c.")
print("From the condition lim_{x->-inf} (f(x) + 3) / e^x = 1:")
print("The numerator must go to 0 as e^x -> 0, so c + 3 = 0, which means c = -3.")
c = -3
print("The limit then simplifies to lim_{x->-inf} (a*e^x + b) = 1, which means b = 1.")
b_coeff = 1

# Now we have f(x) = a*e^(2x) + e^x - 3.
# Use the condition f(ln 2) = 0 to find the coefficient 'a'.
print("From the condition f(ln 2) = 0:")
print("a * e^(2*ln(2)) + e^(ln(2)) - 3 = 0")
print("a * 4 + 2 - 3 = 0")
print("4a - 1 = 0")
a_coeff = 1/4
print(f"So, the coefficient a = {a_coeff}.")
print(f"The function is f(x) = ({a_coeff})e^(2x) + ({b_coeff})e^x + ({c}).")
print("-" * 20)

# Step 2: Analyze the integral equation to find the values of a and b.
print("Step 2: Find the values of a and b from the integral equation.")
print("The integral equation holds if and only if f(ln(b)) = a.")
print(f"Substituting the function f(x), we get the relation:")
print(f"a = ({a_coeff}) * (e^(ln(b)))^2 + ({b_coeff}) * e^(ln(b)) + ({c})")
print(f"a = {a_coeff}*b^2 + {b_coeff}*b {c}")
print("-" * 20)

# Step 3: Solve for a and b.
print("Step 3: To find a unique solution, we test the natural choice a = 0.")
# We set the variable 'a' (the integral limit) to 0.
a_final = 0
print(f"Let a = {a_final}. We solve for b in the equation:")
# 0 = (1/4)b^2 + b - 3
# b^2 + 4b - 12 = 0
coeffs_quadratic = [a_coeff, b_coeff, c]
roots = np.roots(coeffs_quadratic)
print(f"0 = {a_coeff}*b^2 + {b_coeff}*b {c}")
print("This is a quadratic equation. The roots are b = " + " or b = ".join(map(str, roots)))
# The term ln(b) requires b > 0.
b_final = max(roots)
print(f"Since ln(b) is in the equation, b must be positive. So we choose b = {b_final}.")
print(f"The values are a = {a_final} and b = {b_final}.")
print("-" * 20)

# Step 4: Calculate the final answer.
print("Step 4: Calculate the final result.")
result = a_final + b_final
print(f"The value of a + b is {a_final} + {b_final} = {result}.")
