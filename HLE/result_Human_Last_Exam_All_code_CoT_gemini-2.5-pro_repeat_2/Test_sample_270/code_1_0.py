from fractions import Fraction

# This script calculates the fractional Dehn twist coefficient (FDTC)
# based on the algebraic properties of the mapping class group Mod(T_1,1).

# 1. The problem is reduced to an algebraic calculation using a known relation
# in the group: (D_a * D_b)^6 = D_delta^(-1).

# 2. We use the properties of the FDTC function, which we call `f`.
#    - f(g^n) = n * f(g)
#    - f(D_delta) = 1
#    - f applied to the relation gives: 6 * f(D_a * D_b) = -1.

# 3. From the relation, we find the FDTC of the composite twist (D_a * D_b).
f_Da_Db_num = -1
f_Da_Db_den = 6
f_Da_Db = Fraction(f_Da_Db_num, f_Da_Db_den)

# 4. The problem asks for the FDTC of (D_a * D_b)^9.
n = 9
final_result = n * f_Da_Db

# 5. Print the final equation with all its components.
print("The fractional Dehn twist coefficient of the element (D_a * D_b)^9 is calculated as follows:")
print(f"FDTC((D_a * D_b)^9) = {n} * FDTC(D_a * D_b)")
print(f"                      = {n} * ({f_Da_Db_num}/{f_Da_Db_den})")
print(f"                      = {final_result.numerator}/{final_result.denominator}")

print("\nThe final numerical value is:")
print(float(final_result))