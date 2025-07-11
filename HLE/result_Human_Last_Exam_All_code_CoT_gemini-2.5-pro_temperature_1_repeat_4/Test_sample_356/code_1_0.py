import numpy as np

# Step 1: Calculate u1 from the properties of sheets B and C.
# Given u2 = 3. From the center of gravity equations for the welded sheet, we derived: 5/6 * u1 = 8 * u2
u2 = 3.0
u1 = (48.0 / 5.0) * u2

# Step 2: Define the function f(x) and its derivatives f'(x) and f''(x).
def integrand(t):
    """The function inside the integral for f(x)."""
    if 1 + t**4 == 0:
        return 0
    return (2 * t**3 + t) / (1 + t**4)

def first_derivative(x):
    """The first derivative f'(x)."""
    if 1 + x**4 == 0:
        return 0
    return (2 * x**3 + x) / (1 + x**4)

def second_derivative(x):
    """The second derivative f''(x)."""
    denominator = (1 + x**4)**2
    if denominator == 0:
        return 0
    numerator = -2*x**6 - 3*x**4 + 6*x**2 + 1
    return numerator / denominator

# Step 3: Calculate f(5) using Simpson's rule with n=10 subintervals.
n = 10
a_int = 0.0
b_int = 5.0
h = (b_int - a_int) / n
x_vals = np.linspace(a_int, b_int, n + 1)
y_vals = np.array([integrand(val) for val in x_vals])

# Simpson's Rule formula: h/3 * (y0 + 4y1 + 2y2 + ... + 4y(n-1) + yn)
simpson_sum = y_vals[0] + y_vals[-1]
simpson_sum += 4 * np.sum(y_vals[1:-1:2]) # Sum of odd-indexed elements
simpson_sum += 2 * np.sum(y_vals[2:-1:2]) # Sum of even-indexed elements
f5_val = (h / 3.0) * simpson_sum

# Calculate f'(5) and f''(5)
fp5_val = first_derivative(5.0)
fpp5_val = second_derivative(5.0)

# Round the intermediate f-term results to one decimal place as instructed.
f5_rounded = round(f5_val, 1)
fp5_rounded = round(fp5_val, 1)
fpp5_rounded = round(fpp5_val, 1)

# Step 4: Calculate 'a'.
# a = (u1 / 27) * (f(5) - 2*f'(5) + 2*f''(5))^3
expression_val = f5_rounded - 2 * fp5_rounded + 2 * fpp5_rounded
a_val = (u1 / 27.0) * (expression_val**3)

# Step 5: Calculate l.
# From the analysis of sheet A, we found l = 4 * sqrt(3) * a.
l_val = 4 * np.sqrt(3) * a_val

# Final Output Section
print("--- Intermediate Calculations ---")
print(f"Given u2 = {u2}")
print(f"Derived u1 = {u1}")
print(f"f(5) [Simpson's Rule] = {f5_val:.4f}")
print(f"f'(5) = {fp5_val:.4f}")
print(f"f''(5) = {fpp5_val:.4f}")

print("\n--- Rounded f-term values ---")
print(f"f(5)_rounded = {f5_rounded}")
print(f"f'(5)_rounded = {fp5_rounded}")
print(f"f''(5)_rounded = {fpp5_rounded}")

print("\n--- Calculation of 'a' ---")
print(f"a = ({u1} / 27) * ({f5_rounded} - 2 * {fp5_rounded} + 2 * {fpp5_rounded})^3")
print(f"a = {u1 / 27.0:.4f} * {expression_val:.1f}^3")
print(f"Value of a = {a_val}")

print("\n--- Final Calculation of 'l' ---")
print(f"l = 4 * sqrt(3) * a")
print(f"l = {4} * {np.sqrt(3):.4f} * {a_val}")
print(f"Final value of l = {l_val:.4f}")

print(f"\n<<<l = {l_val}>>>")