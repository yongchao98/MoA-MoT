import numpy as np

# From the derivation, any non-zero eigenvalue lambda must satisfy the equation
# lambda^4 = 1. This means the possible non-zero eigenvalues are the 4th roots of unity.
# The problem asks for the number of non-real eigenvalues. We can find these roots by
# solving the polynomial equation z^4 - 1 = 0.
# The coefficients of z^4 + 0*z^3 + 0*z^2 + 0*z - 1 = 0 are [1, 0, 0, 0, -1].
coefficients = [1, 0, 0, 0, -1]
roots = np.roots(coefficients)

# We also found that lambda = 0 is a possible eigenvalue.
# Let's add it to our list of possible eigenvalues.
all_possible_eigenvalues = np.append(roots, 0)

# We are interested in the set S of non-real eigenvalues.
# We filter the list of all possible eigenvalues to find the non-real ones.
non_real_eigenvalues = []
for lam in all_possible_eigenvalues:
    # A complex number is non-real if its imaginary part is non-zero.
    # We use a small tolerance for floating point comparison.
    if abs(lam.imag) > 1e-9:
        non_real_eigenvalues.append(lam)

# The set S consists of these non-real eigenvalues. We want to find its largest possible size.
# This is the number of distinct non-real solutions we found. We can construct a matrix
# that has all of them as eigenvalues.

print("The complete set of possible eigenvalues is: ", {np.round(v, 5) for v in all_possible_eigenvalues})
print("The subset S of non-real eigenvalues is: {", end="")

# The following code is for formatting the output of S nicely
formatted_eigenvalues = []
for val in non_real_eigenvalues:
    if abs(val.real) < 1e-9 and abs(val.imag - 1) < 1e-9:
        formatted_eigenvalues.append("i")
    elif abs(val.real) < 1e-9 and abs(val.imag + 1) < 1e-9:
        formatted_eigenvalues.append("-i")
    else:
        # Fallback for unexpected formats
        formatted_eigenvalues.append(f"{val.real:.2f}{val.imag:+.2f}i")
print(", ".join(formatted_eigenvalues), end="")
print("}")

largest_size_S = len(non_real_eigenvalues)
print(f"The largest possible size of the set S is: |S| = {largest_size_S}")