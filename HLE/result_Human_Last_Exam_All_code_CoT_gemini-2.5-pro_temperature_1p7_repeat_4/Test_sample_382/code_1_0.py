import numpy as np

# 1. Explain the theory and the plan
print("We want to find the greatest possible rank of a matrix E with minimum Frobenius norm")
print("such that a given non-zero vector x is a least-squares solution to min_z ||(A+E)z - b||_2.")
print("\nAccording to matrix perturbation theory, the solution E that minimizes ||E||_F has a rank of at most 1.")
print("A rank of 1 is achievable if the original residual, r = Ax - b, is non-zero.")
print("\nThe perturbation E that makes the new residual zero is given by E = u * v^T, where:")
print("u = -(Ax - b)")
print("v^T = x^T / (x^T * x)")
print("This is an outer product, which results in a rank-1 matrix if u and v are non-zero.")

# 2. Set up and display a concrete example
A = np.array([[1., 2.], [3., 4.], [5., 6.]])
b = np.array([1., 1., 1.])
x = np.array([1., -2.])

print("\nLet's use a specific example:")
print("A = \n", A)
print("b = \n", b)
print("x = \n", x)

# 3. Calculate the components of the formula and output them
r = A @ x - b
u = -r
xTx = x.T @ x
# The formula is E = -outer(r, x) / (x.T @ x). This is outer(u, x) / xTx.
# This is the outer product of u and (x / xTx).
# So v^T is x^T / xTx
v_T = x.T / xTx

print("\nFirst, we calculate r = Ax - b:")
print("r =\n", r)
print("\nThe vector u = -r is:")
print("u =\n", u)
print("\nNext, we calculate x^T * x:")
print(f"x^T * x = {xTx}")
print("\nThe vector v^T = x^T / (x^T * x) is:")
print("v^T =\n", v_T)

# 4. Calculate and display the final matrix E and its rank
E = np.outer(u, v_T)

print("\nThe matrix E is the outer product u * v^T:")
print("E =\n", E)

rank_E = np.linalg.matrix_rank(E)

print(f"\nThe rank of this matrix E is {rank_E}.")
print("\nSince we have shown a case where the rank is 1, and theory confirms the rank is at most 1,")
print("the greatest possible rank is 1.")

<<<1>>>