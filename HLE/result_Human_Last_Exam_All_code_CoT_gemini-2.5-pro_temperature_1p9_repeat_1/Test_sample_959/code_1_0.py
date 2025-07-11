import math

# The problem simplifies to finding the order of the outer automorphism group
# of the cyclic group C of order 31.
# G is the trivial group {1}. Any central extension E of G by C is isomorphic to C.
# So E is isomorphic to the cyclic group of order 31, Z_31.
# The set of extensions E has only one element.
# The sum is the order of the outer automorphism group of C, o(C).
# o(C) = |Out(C)| = |Aut(C) / Inn(C)|.
# Since C is abelian, Inn(C) is trivial.
# So |Out(C)| = |Aut(C)|.
# The order of the automorphism group of a cyclic group of order n is phi(n),
# where phi is Euler's totient function.
# C is of order 31, which is a prime number.
p = 31
# For a prime p, phi(p) = p - 1.
order = p - 1

print(f"{p} - 1 = {order}")
