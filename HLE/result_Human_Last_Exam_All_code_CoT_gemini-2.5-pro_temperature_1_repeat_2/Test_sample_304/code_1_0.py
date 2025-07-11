# Plan:
# 1. The problem is to find the maximum value of H(x, y, z, s1, s2).
# 2. Based on the constraints, we can simplify the expression and find an upper bound.
#    - From H(s1|z,x) = 0 and H(s2|y,z) = 0, we deduce that H(x,y,z,s1,s2) = H(x,y,z).
#    - From H(y|x,s2) = 0 and H(s2|y,z) = 0, we can show that H(y|x,z) = 0.
#    - This simplifies the joint entropy: H(x,y,z) = H(x,z) + H(y|x,z) = H(x,z).
#    - Using the property H(A,B) <= H(A) + H(B) and the constraints H(x)<=1, H(z)<=1,
#      we get H(x,y,z) <= H(x) + H(z) <= 1 + 1 = 2.
# 3. This upper bound of 2 is achievable. A construction is provided in the explanation above.
#    The maximum value is achieved when H(x) and H(z) are maximized and x and z are independent.

# Maximal entropy for x
H_x_max = 1

# Maximal entropy for z
H_z_max = 1

# The maximal joint entropy H(x,y,z,s1,s2) is given by H(x) + H(z) under the condition
# that x and z are independent, which is consistent with the constraints.
max_entropy = H_x_max + H_z_max

print(f"The derivation shows that the maximal entropy is bounded by H(x) + H(z).")
print(f"Given the constraints H(x) <= {H_x_max} and H(z) <= {H_z_max}, the maximal value is achieved when they are independent and take their maximal values.")
print(f"Maximal Entropy = H(x) + H(z) = {H_x_max} + {H_z_max} = {max_entropy}")