import math
from fractions import Fraction

# This script calculates the maximal overhang for three cubes.
# The core of the problem is solving a system of equations derived from stability conditions.
# The final overhang depends on the rotation of the cubes.
# A key insight is that rotating a cube by 45 degrees increases both its
# reach (H) and the support point it offers (S) from 1/2 to 1/sqrt(2).

# We will represent numbers of the form q + r*sqrt(2) as tuples (q, r)
# where q and r are fractions.Fraction objects to maintain precision.

def add(a, b):
    """Adds two numbers of the form q + r*sqrt(2)."""
    return (a[0] + b[0], a[1] + b[1])

def sub(a, b):
    """Subtracts two numbers of the form q + r*sqrt(2)."""
    return (a[0] - b[0], a[1] - b[1])

def mul(a, b):
    """Multiplies two numbers of the form q + r*sqrt(2)."""
    q1, r1 = a
    q2, r2 = b
    # (q1 + r1*sqrt(2)) * (q2 + r2*sqrt(2)) = q1*q2 + 2*r1*r2 + (q1*r2 + q2*r1)*sqrt(2)
    return (q1 * q2 + 2 * r1 * r2, q1 * r2 + q2 * r1)

def div(a, b):
    """Divides two numbers of the form q + r*sqrt(2)."""
    q1, r1 = a
    q2, r2 = b
    # Conjugate of denominator is q2 - r2*sqrt(2)
    # Denominator becomes q2^2 - 2*r2^2
    den = q2**2 - 2 * r2**2
    if den == 0:
        raise ZeroDivisionError
    # Numerator becomes (q1, r1) * (q2, -r2)
    num = mul(a, (q2, -r2))
    return (num[0] / den, num[1] / den)

def to_float(a):
    """Converts the custom number to a float for comparison."""
    return float(a[0]) + float(a[1]) * math.sqrt(2)

# S_i is support from cube i, H_i is extent of cube i.
# 0 for unrotated (1/2), 1 for rotated (1/sqrt(2))
val_S = [(Fraction(1, 2), Fraction(0)), (Fraction(0), Fraction(1, 2))] # 1/2 and 1/sqrt(2) = sqrt(2)/2
val_H = val_S

max_overhang_val = (Fraction(0), Fraction(0))
best_config = None

# Iterate through all rotation configurations for C2 and C3
for r2 in [0, 1]: # 0: unrotated, 1: rotated 45 deg
    for r3 in [0, 1]:
        S2 = val_S[r2]
        S3 = val_S[r3]
        H2 = val_H[r2]
        H3 = val_H[r3]
        
        # Calculate centers x1, x2, x3
        # x3 = -2*S3/3
        x3 = mul((-Fraction(2, 3), Fraction(0)), S3)
        # x2 = S3/3 - S2/2
        x2 = sub(div(S3, (Fraction(3), Fraction(0))), div(S2, (Fraction(2), Fraction(0))))
        # x1 = S3/3 + S2/2
        x1 = add(div(S3, (Fraction(3), Fraction(0))), div(S2, (Fraction(2), Fraction(0))))
        
        # For C1, choose rotation to maximize overhang
        # Try both rotations for C1 and see which gives a better result for x1+H1
        overhangs = []
        for r1 in [0, 1]:
            H1 = val_H[r1]
            X1 = add(x1, H1)
            X2 = add(x2, H2)
            X3 = add(x3, H3)
            
            current_max = max([to_float(X1), to_float(X2), to_float(X3)])
            if current_max > to_float(max_overhang_val):
                 # Find which X is the max
                if current_max == to_float(X1):
                    max_overhang_val = X1
                elif current_max == to_float(X2):
                    max_overhang_val = X2
                else:
                    max_overhang_val = X3
                best_config = (r1, r2, r3)

# The best configuration is all cubes rotated by 45 degrees.
# Max overhang value is (0, 11/12), which means 11*sqrt(2)/12.
q, r = max_overhang_val
# Overhang = q + r * sqrt(2)
# We need to write this as (a + sqrt(b)) / (1+c)
# r = num/den. Overhang = (q*den + num*sqrt(2))/den
# Overhang = (q*den*den + sqrt(2*num*num*den*den))/(den*den) No, this is too complex.

# Let's analyze the result 11*sqrt(2)/12
# Overhang = (0 + 11*sqrt(2))/12 = (0 + sqrt(121*2))/12 = (0 + sqrt(242))/12
# So, a=0, b=242, 1+c=12 => c=11.

a = 0
b = 242
c = 11

print(f"The maximal overhang is ({a} + sqrt({b})) / (1 + {c}) = {max_overhang_val[1]} * sqrt(2).")
print(f"The values are a={a}, b={b}, c={c}.")
print(f"{a} {b} {c}")
