# This script requires a SageMath environment to run.
# You can use it in a SageMath notebook or by installing SageMath as a Python package.

from sage.all import EllipticCurve, DirichletGroup, QQbar, factorial, round

# 1. Define the Elliptic Curve
# The curve is given by y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
# For y^2 + y = x^3 - x^2 - 10x - 20, the coefficients are:
# a1=0, a2=-1, a3=1, a4=-10, a6=-20
E = EllipticCurve([0, -1, 1, -10, -20])

# 2. Compute the rank r of the Mordell-Weil group E(Q)
r = E.rank()

# 3. Find the primitive cubic Dirichlet characters of conductor 7
# We use QQbar as the codomain to get complex values for the characters.
G = DirichletGroup(7, codomain=QQbar)
cubic_chars = [chi for chi in G if chi.order() == 3 and chi.is_primitive()]
chi1 = cubic_chars[0]
chi2 = cubic_chars[1]

# 4. Compute the leading coefficients a and b
# Get the L-series for each twisted curve
L1 = E.lseries().twist(chi1)
L2 = E.lseries().twist(chi2)

# Determine the order of vanishing (analytic rank) at s=1
rank1 = L1.rank()
rank2 = L2.rank()

# The leading coefficient of the Taylor series is L^(k)(1) / k! where k is the analytic rank.
# Sage's central_value() function returns L^(k)(1).
a = L1.central_value() / factorial(rank1)
b = L2.central_value() / factorial(rank2)

# 5. Calculate the final sum and print the results
# a and b are complex conjugates, so their sum is a real number.
total_sum = r + a + b

# Print the components of the final equation
# We format the complex numbers a and b for clear output.
# The .n(digits=10) provides a numerical approximation.
a_str = f"({a.real().n(digits=8)} + {a.imag().n(digits=8)}i)"
b_str = f"({b.real().n(digits=8)} + {b.imag().n(digits=8)}i)"

print(f"The rank of E(Q) is r = {r}")
print(f"The leading coefficient for L(E, s, chi1) is a = {a_str}")
print(f"The leading coefficient for L(E, s, chi2) is b = {b_str}")

# The final result rounded to four decimal places
final_answer = round(total_sum, 4)

# Print the final equation and the result
print(f"\nThe final equation is r + a + b = {r} + {a_str} + {b_str}")
print(f"The sum r + a + b is approximately {total_sum.n(digits=10)}")
print(f"The result rounded to four decimal places is: {final_answer}")