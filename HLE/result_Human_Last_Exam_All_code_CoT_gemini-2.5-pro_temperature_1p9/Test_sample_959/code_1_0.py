import math

# The problem is to compute the sum of orders of outer automorphism groups
# for all central extensions E of a group G by C.
#
# Our analysis concluded that the group G is trivial.
# G = {1}.
# The set of central extensions E has only one element, which is C itself.
# E = C.
# So the sum is just the order of the outer automorphism group of C.
#
# Let C be the cyclic group of order p.
p = 31

# For an abelian group, the outer automorphism group is the same as the automorphism group.
# The order of the automorphism group of a cyclic group of order p (a prime) is phi(p).
# phi(p) = p - 1 for a prime p.
order_aut_C = p - 1

# The final sum is this single value.
sum_of_orders = order_aut_C

print(f"Let C be the cyclic group of order p.")
print(f"p = {p}")
print(f"The order of the outer automorphism group of C is o(C) = p - 1.")
print(f"o(C) = {p} - 1 = {order_aut_C}")
print(f"The sum is therefore equal to {sum_of_orders}.")
