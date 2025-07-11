import cmath

# The final analytical expression for the integral is given in terms of the Hurwitz zeta function.
# Let I be the value of the integral.
# I = i * sqrt(pi) * [zeta(3/2, 3/2 - i*2/3) + zeta(3/2, 3/2 + i*2/3)]
# This expression represents the analytical value. We can print this relationship.

# Parameters for the final expression
s = 3/2
a = 3/2
b = 2/3

# We will print the equation that gives the value of the integral.
# Note that we format the numbers to match the derivation.
print("The analytical value of the integral I is given by the equation:")
print(f"I = i * sqrt(pi) * [zeta({s}, {a} - {b}*i) + zeta({s}, {a} + {b}*i)]")

# This formula can be numerically evaluated to high precision.
# For example, using the mpmath library:
from mpmath import mp, j, pi, sqrt, zeta
mp.dps = 30 # set precision
s_mp = mp.mpf('1.5')
a_mp = mp.mpf('1.5')
b_mp = mp.mpf('2')/mp.mpf('3')

val = j * sqrt(pi) * (zeta(s_mp, a_mp - j*b_mp) + zeta(s_mp, a_mp + j*b_mp))
print("\nA numerical evaluation of this expression yields:")
print(f"I ≈ {val}")
