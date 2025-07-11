import math

def is_perfect_square(n):
    """Checks if a number is a positive perfect square."""
    if n <= 0:
        return False
    sqrt_n = int(math.sqrt(n))
    return sqrt_n * sqrt_n == n

def count_parallelograms():
    """
    Counts the number of distinct parallelograms based on the problem's restrictions.
    """
    found_parallelograms = set()

    # 1. Iterate through side lengths a and b.
    # From 2a < a+b < 100, we get a < b and a+b < 100.
    for a in range(1, 50):
        for b in range(a + 1, 100 - a):
            
            # 2. Check restrictions on a and b.
            # Must be coprime.
            if math.gcd(a, b) != 1:
                continue
            
            # For diagonals to be integers, a and b must have opposite parity.
            if a % 2 == b % 2:
                continue

            ab = a * b
            sum_of_squares = a**2 + b**2

            # 3. Find all Pythagorean triples with hypotenuse 'ab'.
            # Generate primitive triples (x, y, z) and scale them.
            # z = m^2 + n^2. Limit m to where m^2 is not much larger than ab.
            limit_m = int(math.sqrt(ab)) + 1
            for m in range(2, limit_m):
                for n in range(1, m):
                    # For primitive triples, m and n are coprime with opposite parity.
                    if (m - n) % 2 == 1 and math.gcd(m, n) == 1:
                        z = m**2 + n**2
                        
                        if ab % z == 0:
                            # We found a basis for a valid triple.
                            k = ab // z
                            # Legs of the scaled triple.
                            leg1 = k * (m**2 - n**2)
                            leg2 = k * (2 * m * n)

                            # 'Y' can be either leg.
                            for Y in {leg1, leg2}:
                                if Y == 0:  # Excludes rectangles
                                    continue
                                
                                # 4. Check if diagonals are integers.
                                d1_sq = sum_of_squares - 2 * Y
                                d2_sq = sum_of_squares + 2 * Y

                                if is_perfect_square(d1_sq) and is_perfect_square(d2_sq):
                                    # Found a valid parallelogram, defined by (a, b, Y).
                                    found_parallelograms.add((a, b, Y))

    return len(found_parallelograms)

# Execute the search and print the final count.
result = count_parallelograms()
print(result)