import cmath

# Based on the mathematical derivation, any eigenvalue lambda of a matrix A
# satisfying A^3 = A^* must solve the equation: lambda^3 = conjugate(lambda).
# The full set of solutions to this equation is {0, 1, -1, i, -i}.

all_possible_eigenvalues = {0, 1, -1, 1j, -1j}

# The problem asks for the largest size |S| of a set S of non-real eigenvalues.
# These are complex numbers with a non-zero imaginary part.
# We can filter the set of all possible eigenvalues to find them.
S = {val for val in all_possible_eigenvalues if val.imag != 0}

print("The set S consists of the non-real eigenvalues that satisfy the equation lambda^3 = conjugate(lambda).")
print("The elements of this set are:")

# The final equation is S = {i, -i}. We print each number in this set.
for number in S:
    print(number)

size_of_S = len(S)
print(f"\nThe largest size of such a set S is {size_of_S}.")