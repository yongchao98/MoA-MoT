import math

def is_perfect_square(n):
    """Checks if a non-negative integer is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def solve():
    """
    Calculates the number of distinct parallelograms based on the problem's restrictions.
    """
    found_parallelograms = set()
    limit_sum = 100

    # Iterate through possible side lengths a and b
    # Condition 2: 2a < a + b < 100 implies a < b and a + b < 100
    for a in range(1, limit_sum // 2):
        for b in range(a + 1, limit_sum - a):
            # Condition 2: a and b are coprime
            if math.gcd(a, b) != 1:
                continue

            # Condition 4: Diagonals are integers. From the parallelogram law, d1^2 + d2^2 = 2(a^2 + b^2).
            N = 2 * (a*a + b*b)

            # Find integer pairs (d1, d2) that solve the equation.
            # We assume d1 < d2 to satisfy Condition 1 (not a rectangle).
            # This implies d1^2 < N/2.
            d1_limit = int(math.sqrt(N / 2))
            for d1 in range(1, d1_limit + 1):
                d2_squared = N - d1*d1
                if is_perfect_square(d2_squared):
                    d2 = int(math.sqrt(d2_squared))

                    # A non-degenerate parallelogram requires its constituent triangles to be non-degenerate.
                    # The triangle inequalities simplify to d1 > b-a and a+b > d2.
                    if not (d1 > b - a and a + b > d2):
                        continue

                    # Condition 3: The area of the parallelogram must be an integer.
                    # Area = sqrt(P)/2, where P = (d1-(b-a))(d1+(b-a))(a+b-d1)(a+b+d1).
                    # For the area to be an integer, P must be a perfect square, and its square root must be even.
                    P = (d1 - (b - a)) * (d1 + (b - a)) * ((a + b) - d1) * ((a + b) + d1)

                    if P > 0 and is_perfect_square(P):
                        k = int(math.sqrt(P))
                        if k % 2 == 0:
                            # All conditions met. We have found a distinct parallelogram.
                            # We add the canonical tuple (a, b, d1, d2) to a set to count unique solutions.
                            found_parallelograms.add((a, b, d1, d2))

    # The total number of distinct parallelograms is the size of the set.
    print(len(found_parallelograms))

solve()