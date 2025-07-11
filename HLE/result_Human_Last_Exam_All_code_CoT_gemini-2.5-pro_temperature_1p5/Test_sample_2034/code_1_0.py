# The problem is to determine the coefficients c1 and c2 in the Jaakola-Haussler bound
# extension for beta-SVMs.
# The bound is of the form:
# -(K α^{D-i})_i ≤ +(1 + c1*β)α^{D}_i - (1 + c2*β)(K α^{D})_i + o(β)
#
# Based on the derivation, we found the values for c1 and c2.
# Here we will print these values and the numbers that appear in the final equation.

c1 = -1
c2 = 0

print("The determined coefficients are:")
print(f"c1 = {c1}")
print(f"c2 = {c2}")

print("\nThe final extended bound is:")
print(f"-(K α^D-i)_i ≤ (1 + ({c1})β)α^D_i - (1 + ({c2})β)(K α^D)_i + o(β)")
print(f"-(K α^D-i)_i ≤ (1 - β)α^D_i - (K α^D)_i + o(β)")


print("\nThe numbers in the final equation +(1 + c1*β)α - (1 + c2*β)(Kα) are:")
# The numbers are 1, c1, 1, and c2.
print(1)
print(c1)
print(1)
print(c2)