import math

# The number of Fq-rational maximal tori of a group of type E8 is equal to the
# number of conjugacy classes in its Weyl group, W(E8).

# For context, we first compute the order of W(E8).
# This is calculated from the degrees of the fundamental invariants of the group.
# The degrees (d_i) are related to the exponents (m_i) by the formula d_i = m_i + 1.

# The exponents for the E8 root system are:
exponents_E8 = [1, 7, 11, 13, 17, 19, 23, 29]

# The degrees of the fundamental invariants are:
degrees_E8 = [m + 1 for m in exponents_E8]

# The order of W(E8) is the product of these degrees.
order_W_E8 = math.prod(degrees_E8)

print("For context, let's compute the order of the Weyl group W(E8).")
print(f"The degrees of the fundamental invariants of E8 are: {degrees_E8}")
degree_str = " * ".join(map(str, degrees_E8))
print(f"The order |W(E8)| is the product of the degrees:")
print(f"|W(E8)| = {degree_str} = {order_W_E8}")
print("-" * 30)

# The number of conjugacy classes of W(E8) gives the number of maximal tori.
# This result is a known fact from the theory of reflection groups.
num_tori = 112

# The final answer is given by the following equation:
print("Number of F_q-rational maximal tori in G(E8) = Number of conjugacy classes in W(E8)")
print(f"Number of conjugacy classes in W(E8) = {num_tori}")
