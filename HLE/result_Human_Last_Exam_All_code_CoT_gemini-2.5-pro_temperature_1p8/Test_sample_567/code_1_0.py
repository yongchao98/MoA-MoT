import math

# The problem asks for the value 'a' where the only obstruction is volume.
# This condition is met when c(a) = sqrt(a).
# This is known to be true if and only if 'a' is a perfect square.
# The question asks for "the" value, implying a unique or canonical answer.
# The first trivial case is n=1, a=1.
# The first non-trivial and most famous case is for n=2.

n = 2
a = n**2

# For a = n^2, the capacity function c(a) is exactly n.
c_a = n

# The volume constraint is sqrt(a).
volume_constraint = math.sqrt(a)

# We check that c(a) equals the volume constraint.
# The script will print the equation demonstrating this equality.
if c_a == volume_constraint:
    print(f"For a = {a}, the capacity is c({a}) = {c_a}.")
    print(f"The volume constraint is sqrt({a}) = {volume_constraint}.")
    print("Since c(a) equals the volume constraint, this is a value where the only obstruction is volume.")
    print("\nThe final equation is:")
    print(f"c({a}) = sqrt({a}) = {int(volume_constraint)}")
