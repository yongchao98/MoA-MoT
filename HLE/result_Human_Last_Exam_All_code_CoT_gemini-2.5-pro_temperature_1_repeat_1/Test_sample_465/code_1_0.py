import numpy as np
import cmath

# The chromatic polynomial of the 4-cycle graph (C4) is P(q) = q^4 - 4q^3 + 6q^2 - 3q.
# We want to find the roots of the equation P(q) = 0.
# The coefficients of the polynomial, from the highest power to the lowest, are [1, -4, 6, -3, 0].
coeffs = [1, -4, 6, -3, 0]

# Let's print the equation we are solving.
# The prompt asked to output each number in the final equation.
print("The chromatic polynomial equation for C4 is:")
print(f"({coeffs[0]})q^4 + ({coeffs[1]})q^3 + ({coeffs[2]})q^2 + ({coeffs[3]})q + ({coeffs[4]}) = 0")
print("\nFinding the roots of this polynomial:")

# Use numpy's `roots` function to find the roots of the polynomial.
roots = np.roots(coeffs)

# Sort the roots for cleaner presentation
# Custom sort key to handle complex numbers
def sort_key(x):
    return (x.real, x.imag)
sorted_roots = sorted(roots, key=sort_key)


print("The chromatic roots are:")
for root in sorted_roots:
    # Format complex numbers nicely
    if abs(root.imag) < 1e-9:
        print(f"{root.real:.4f}")
    else:
        print(f"{root.real:.4f} + {root.imag:.4f}i")

print("\nAs we can see, two of the roots are complex (non-real), which confirms statement B.")
