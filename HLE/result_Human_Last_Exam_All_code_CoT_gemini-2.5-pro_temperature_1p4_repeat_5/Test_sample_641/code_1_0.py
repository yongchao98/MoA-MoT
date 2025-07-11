# The prime field order is q=997. The group is PSU(4, 997).
q = 997

# N1 is the number of non-central involutions in SU(4, q).
# An involution is an element g != I such that g^2 = I.
# The formula for this number is q^4 * (q^2 - q + 1) * (q^2 + 1).
n1 = q**4 * (q**2 - q + 1) * (q**2 + 1)

# N2 is the number of quasi-involutions in SU(4, q).
# A quasi-involution is an element g such that g^2 = -I.
# For n=4 and q=1(mod 4), the formula is q^5 * (q^3 + 1) * (q^4 - 1).
n2 = q**5 * (q**3 + 1) * (q**4 - 1)

# The total number of involutions in PSU(4, q) is the sum of involutions
# coming from both types of elements, divided by the order of the center (which is 2).
# This gives the final formula: (N1 + N2) / 2.
total_involutions = (n1 + n2) // 2

# The problem asks to output each number in the final equation.
# The final equation is (n1 + n2) / 2 = total_involutions.
print(f"({n1} + {n2}) // 2 = {total_involutions}")
