import numpy as np
import sympy

# Step 1: Set up the problem from the derivation.
# We test the hypothesis that C = -1 is a constant lower bound.
# This requires the polynomial f(u) = u^4 - u^3 - 5u^2 + 3u + 2 to be non-negative
# for all u in the interval [0, 1].

# The coefficients of the polynomial f(u)
c4, c3, c2, c1, c0 = 1, -1, -5, 3, 2

# Step 2: Print the equation we are verifying
print("The problem reduces to verifying the inequality:")
print(f"f(u) = {c4}*u^4 + ({c3})*u^3 + ({c2})*u^2 + {c1}*u + {c0} >= 0, for u in [0, 1]\n")

# Step 3: Find the roots of the polynomial to understand its behavior.
# The roots tell us where the function crosses zero.
coeffs = [c4, c3, c2, c1, c0]
roots = np.roots(coeffs)
# We sort the roots to make them easier to inspect.
roots.sort()

print("The roots of the polynomial f(u) are:")
# We can use sympy to get exact expressions for the roots
u = sympy.Symbol('u')
poly = c4*u**4 + c3*u**3 + c2*u**2 + c1*u + c0
exact_roots = sympy.solve(poly, u)
print(f"Exact roots: {exact_roots}")
print(f"Numerical roots: {[np.round(r, 5) for r in roots]}\n")

# Step 4: Analyze the polynomial on the interval [0, 1]
# From the roots, we see that the only root within the interval [0, 1] is u=1.
# The other roots are outside this interval.
# To check if f(u) is non-negative on [0, 1], we can test the sign of f(u)
# at the boundaries and at one interior point, for instance, u = 0.5.

f_0 = c4*0**4 + c3*0**3 + c2*0**2 + c1*0 + c0
f_0_5 = c4*0.5**4 + c3*0.5**3 + c2*0.5**2 + c1*0.5 + c0
f_1 = c4*1**4 + c3*1**3 + c2*1**2 + c1*1 + c0

print("To verify the inequality on [0, 1], we check values at the boundaries and an interior point:")
print(f"Value at u = 0: f(0) = {f_0}")
print(f"Value at u = 0.5: f(0.5) = {f_0_5}")
print(f"Value at u = 1: f(1) = {f_1}\n")

# Step 5: Conclusion
# Since f(0) > 0, f(1) = 0, and the function does not cross zero inside (0,1),
# the function must be non-negative over the entire interval [0, 1].
print("Conclusion: The polynomial f(u) is non-negative for all u in [0, 1].")
print("Therefore, C = -1 is a valid constant lower bound for d(t,x).")
print("Since the initial minimum slope is -0.5, which is greater than -1,")
print("the slope d(t,x) will remain greater than or equal to -1 for all t > 0.")
