import math

def solve():
    """
    This function finds the number of ordered triplets of positive integers (x, y, z)
    such that x*y*z = 216 and x^(log2(y)) = y^(log2(z)).
    
    The second condition simplifies to (y=1 or x=z). We use the Principle of
    Inclusion-Exclusion to find the total number of solutions.
    """
    volume = 216

    # Let A be the set of solutions where y = 1.
    # The equation becomes x * z = 216. The number of solutions equals the number of divisors of 216.
    count_A = 0
    for i in range(1, volume + 1):
        if volume % i == 0:
            count_A += 1

    # Let B be the set of solutions where x = z.
    # The equation becomes x^2 * y = 216. We need to find the number of possible integer values for x.
    # This means x^2 must be a divisor of 216.
    count_B = 0
    limit = math.isqrt(volume)
    for x in range(1, limit + 1):
        # We check if x*x is a divisor of the volume.
        if volume % (x * x) == 0:
            count_B += 1

    # Let A_intersect_B be the set of solutions where y=1 AND x=z.
    # The equation becomes x^2 = 216.
    count_AB = 0
    # We check if volume is a perfect square.
    if limit * limit == volume:
        count_AB = 1
    
    # The total number of solutions is N(A) + N(B) - N(A and B).
    total_permutations = count_A + count_B - count_AB
    
    print(f"Number of permutations where y = 1: {count_A}")
    print(f"Number of permutations where x = z: {count_B}")
    print(f"Number of permutations where y = 1 AND x = z: {count_AB}")
    print(f"Total number of permutations = {count_A} + {count_B} - {count_AB} = {total_permutations}")

solve()