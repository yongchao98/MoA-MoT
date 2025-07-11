import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    x = int(math.sqrt(n))
    return x * x == n, x

def solve():
    """
    Finds the number of distinct parallelograms based on the given restrictions.
    """
    # A set to store unique parallelograms found, represented by a tuple
    # of their side and diagonal lengths (a, b, m, n) with a < b and m < n.
    found_parallelograms = set()

    # 1. Iterate through all possible side lengths a and b.
    # The condition a + b < 100 with a < b implies a can go up to 49.
    for a in range(1, 50):
        # b must be greater than a and a + b < 100.
        for b in range(a + 1, 100 - a):
            
            # 2. Check side length restrictions.
            # a != b is guaranteed by the loop ranges.
            # 2a < a + b < 100 is also guaranteed by the loop ranges.
            if gcd(a, b) != 1:
                continue

            # 3. Find integer diagonals m and n.
            # From the parallelogram law, m^2 + n^2 = 2(a^2 + b^2).
            K = 2 * (a**2 + b**2)

            # To find solutions for m^2 + n^2 = K, we iterate m up to sqrt(K/2).
            # This ensures m <= n.
            m_limit = int(math.sqrt(K / 2))
            for m in range(1, m_limit + 1):
                n_squared = K - m**2
                is_sq, n = is_perfect_square(n_squared)

                if is_sq:
                    # Found an integer diagonal pair (m, n).
                    
                    # 4. Check remaining parallelogram restrictions.
                    # Restriction: Must not be a rectangle (m != n).
                    # The loop condition m <= m_limit ensures m <= n.
                    # If m == n, we skip.
                    if m == n:
                        continue

                    # Restriction: Diagonals must form valid triangles with the sides.
                    # This means b-a < m < a+b and b-a < n < a+b.
                    # Since m < n, we only need to check m against the lower bound
                    # and n against the upper bound.
                    if not (b - a < m and n < a + b):
                        continue
                        
                    # Restriction: Area must be an integer.
                    # This is true if ((a+b)^2 - m^2) * (m^2 - (b-a)^2) is a perfect square.
                    P = ((a + b)**2 - m**2) * (m**2 - (b - a)**2)
                    
                    # The triangle inequality check above ensures P > 0.
                    is_area_sq, _ = is_perfect_square(P)

                    if is_area_sq:
                        # All conditions met. This is a valid parallelogram.
                        parallelogram = (a, b, m, n)
                        if parallelogram not in found_parallelograms:
                            found_parallelograms.add(parallelogram)
                            print(f"Found parallelogram with sides (a={a}, b={b}) and diagonals (d1={m}, d2={n}).")

    count = len(found_parallelograms)
    print(f"\nTotal number of distinct parallelograms found: {count}")
    return count

# Execute the solver function.
solve()