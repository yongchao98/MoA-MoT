# The problem reduces to computing the order of the automorphism group of C_31.
# For a cyclic group of prime order p, the order of its automorphism group is phi(p) = p - 1.
p = 31
order = p - 1

# The collection of central extensions has only one element, E.
# The order of the outer automorphism group of E is o(E).
# The sum is over a single element, so the sum is just o(E).
# We found o(E) = 30.
# The final equation is simply the calculation of phi(31).
print(f"The order of the automorphism group is phi({p}) = {p} - 1 = {order}")
print(f"The sum is {order}")
