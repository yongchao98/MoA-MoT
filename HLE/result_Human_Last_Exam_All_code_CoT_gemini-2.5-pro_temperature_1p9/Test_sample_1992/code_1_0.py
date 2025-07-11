import numpy as np

# The problem of finding Quantum Temporal Fixed Points (QTFPs) reduces to
# solving a polynomial equation for the proposition's probability 'p'.
# The derivation results in the equation: 4*p^2 - 4*p + 1 = 0.

# We define the coefficients of the quadratic equation a*p^2 + b*p + c = 0
a = 4
b = -4
c = 1
coefficients = [a, b, c]

# We solve for the roots of the polynomial.
roots = np.roots(coefficients)

# A proposition's probability 'p' must be a real number between 0 and 1.
# We count how many of the solutions meet this criterion.
# np.isreal checks if the number is real (non-complex).
valid_qtfp_count = 0
solution_value = 0

print(f"The equation for the truth probability 'p' of a QTFP is: {a}p^2 + {b}p + {c} = 0")

for root in roots:
    # A root must be a real number and its value must be in the [0, 1] range
    # for it to be a valid probability.
    if np.isreal(root) and 0 <= root <= 1:
        valid_qtfp_count += 1
        # In this specific case, we know there's only one root.
        solution_value = np.real(root)

print(f"The solution for 'p' is: {solution_value}")
print(f"The number of valid solutions for p, and thus the number of quantum temporal fixed points, is:")
print(valid_qtfp_count)

<<<1>>>