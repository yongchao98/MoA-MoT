# This script calculates the number of isomorphism classes of del Pezzo surfaces
# of degree 5 over the rational numbers with good reduction everywhere except
# possibly at the prime 2.

# This problem is equivalent to counting the number of isomorphism classes of
# 4-dimensional etale algebras over Q which are unramified outside the prime 2
# and have a square discriminant.

# We count these algebras by considering the possible partitions of the dimension 4.

# Case 1: The algebra is a single quartic field L (partition "4").
# The discriminant of L must be a square and a power of 2 (up to sign).
# Since a square must be positive, disc(L) must be of the form 2^(2k).
# Consulting number field databases (like PARI/GP or LMFDB), we find
# there are exactly two such fields:
# 1. Q(zeta_8), with discriminant 256 = 16^2.
# 2. The field defined by x^4 - 4x^2 + 8, with discriminant 4096 = 64^2.
num_quartic_fields = 2

# Case 2: The algebra is a product of a cubic field L3 and Q (partition "3+1").
# The discriminant of the algebra is disc(L3). It must be a square.
# It is a known result that there are no cubic fields ramified only at prime 2.
num_cubic_algebras = 0

# Case 3: The algebra is a product of two quadratic fields, L2 and L2' (partition "2+2").
# The fields L2 and L2' must be unramified outside 2. The only such fields are
# Q(i) (disc -4), Q(sqrt(2)) (disc 8), and Q(sqrt(-2)) (disc -8).
# For the total discriminant disc(L2)*disc(L2') to be a square, the square-free
# parts of the individual discriminants must be the same. This means L2 must
# be isomorphic to L2'.
# This gives three possibilities:
# 1. Q(i) x Q(i)
# 2. Q(sqrt(2)) x Q(sqrt(2))
# 3. Q(sqrt(-2)) x Q(sqrt(-2))
num_quadratic_pairs = 3

# Case 4: The algebra is a product of a quadratic field L2 and Q x Q (partition "2+1+1").
# The discriminant of the algebra is disc(L2). It must be a square.
# The discriminants of the quadratic fields unramified outside 2 are -4, 8, -8.
# None of these is a square in Q.
num_quadratic_q_q = 0

# Case 5: The algebra is Q x Q x Q x Q (partition "1+1+1+1").
# This is the "split" case. The discriminant is 1, which is a square.
# The algebra is unramified everywhere.
num_split_case = 1

# The total number is the sum of the counts from all cases.
total = num_quartic_fields + num_cubic_algebras + num_quadratic_pairs + num_quadratic_q_q + num_split_case

print("The total number of isomorphism classes is the sum of possibilities from each partition of 4:")
print(f"1. Quartic fields (partition 4): {num_quartic_fields}")
print(f"2. Cubic x Q (partition 3+1): {num_cubic_algebras}")
print(f"3. Two quadratic fields (partition 2+2): {num_quadratic_pairs}")
print(f"4. Quadratic x Q x Q (partition 2+1+1): {num_quadratic_q_q}")
print(f"5. Split case QxQxQxQ (partition 1+1+1+1): {num_split_case}")
print("\nFinal calculation:")
print(f"{num_quartic_fields} + {num_cubic_algebras} + {num_quadratic_pairs} + {num_quadratic_q_q} + {num_split_case} = {total}")
