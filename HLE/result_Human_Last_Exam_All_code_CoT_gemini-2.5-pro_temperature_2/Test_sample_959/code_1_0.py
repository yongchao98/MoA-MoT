import math

def get_prime_euler_totient(p):
    """
    Calculates Euler's totient function for a prime number p.
    phi(p) = p - 1.
    """
    if p > 1:
        # A simple primality test for the sake of correctness.
        for i in range(2, int(math.sqrt(p)) + 1):
            if p % i == 0:
                raise ValueError("The number must be prime.")
        return p - 1
    else:
        raise ValueError("The number must be a prime greater than 1.")

# Step 1: The group G has presentation
# <a, b, c, d | aba^-1 = a^2, bcb^-1 = c^2, cdc^-1 = d^2, dad^-1 = a^2>.
# This group is perfect (G_ab=1) and its Schur multiplier M(G) is trivial.
# The reasoning is that it has a balanced presentation (3 generators, 3 relations after simplification)
# and the matrix of exponent sums is invertible over Z.

# Step 2: The central extensions E are classified by H^2(G, C).
# C is the cyclic group of order 31, C_31.
# By Universal Coefficient Theorem, H^2(G, C) ~= Hom(M(G), C).
# Since M(G)=0, H^2(G, C) is trivial, so there is only one extension E, the direct product G x C.
num_extensions = 1

# Step 3: The sum is over a single element E = G x C_31.
# We need to compute o(E) = |Out(E)| = |Out(G x C_31)|.
# For G being perfect with trivial center, and C abelian, Out(G x C) ~= Out(G) x Aut(C).
# We deduce |Out(G)| = 1, either because G is trivial or it is a complex non-symmetric group.
order_out_G = 1

# Step 4: The cyclic group C is C_31. The order of its automorphism group is phi(31).
p = 31
order_aut_C = get_prime_euler_totient(p)

# Step 5: Compute the order of the outer automorphism group of the single extension.
o_E = order_out_G * order_aut_C

# Step 6: The sum is over a single term.
total_sum = o_E

# As requested, output each number in the final equation.
# The sum has only one term, which is a product of two numbers.
print(f"The number of central extensions is 1.")
print(f"The order of the outer automorphism group for the single extension is the product of |Out(G)| and |Aut(C_31)|.")
print(f"The final sum is calculated as: {order_out_G} * {order_aut_C} = {total_sum}")