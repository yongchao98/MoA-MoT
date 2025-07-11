import math

def get_gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return abs(a)

def is_perfect_square(n):
    """Checks if a number is a perfect square and returns (bool, sqrt_int)."""
    if n < 0:
        return False, 0
    if n == 0:
        return True, 0
    x = int(math.sqrt(n))
    return x * x == n, x

def solve():
    """
    Finds the number of distinct parallelograms satisfying the given criteria.
    """
    count = 0
    # A set to store unique parallelograms, defined by (a, b, d1, d2)
    # where a < b and d1 < d2 to ensure uniqueness.
    found_parallelograms = set()

    # Iterate through side lengths a and b
    # Constraint: 2a < a + b < 100
    # This simplifies to a < b and a + b < 100
    for a in range(1, 50):
        # To satisfy a < b, b starts from a + 1
        # To satisfy a + b < 100, b is at most 99 - a
        for b in range(a + 1, 100 - a):

            # 1. Sides a and b must be coprime
            if get_gcd(a, b) != 1:
                continue

            # 2. Find potential integer diagonals d1, d2
            # d1^2 + d2^2 = 2*(a^2 + b^2)
            sum_of_sq_diagonals = 2 * (a**2 + b**2)

            # Find pairs (d1, d2) that satisfy the sum of squares
            # To iterate efficiently, d1 <= sqrt(sum_of_sq_diagonals / 2)
            limit = int(math.sqrt(sum_of_sq_diagonals / 2.0))
            for d1 in range(1, limit + 1):
                d2_sq = sum_of_sq_diagonals - d1**2
                is_sq, d2 = is_perfect_square(d2_sq)

                if is_sq and d1 < d2: # d1 < d2 ensures it's not a rectangle
                    
                    # 3. Check for integer area
                    # k^2 + (2*Area)^2 = (2ab)^2
                    # k = (d2^2 - d1^2)/2
                    k_num = d2**2 - d1**2
                    if k_num % 2 != 0:
                        continue # k must be an integer
                    k = k_num // 2
                    
                    two_ab_sq = (2 * a * b)**2
                    if k**2 >= two_ab_sq:
                        continue
                    
                    two_area_sq = two_ab_sq - k**2
                    is_area_sq, two_area = is_perfect_square(two_area_sq)

                    # 2*Area must be even, so Area is integer.
                    if is_area_sq and two_area % 2 == 0:
                        area = two_area // 2
                        
                        # Store the canonical representation of the parallelogram
                        # and check for uniqueness
                        current_parallelogram = (a, b, d1, d2)
                        if current_parallelogram not in found_parallelograms:
                            found_parallelograms.add(current_parallelogram)
                            # The problem asks to output numbers from the "final equation"
                            # This is interpreted as printing the parameters of the found parallelograms.
                            # For example, side^2 + side^2 = diag^2/2 + diag^2/2 or Area formula.
                            # We print the components: a, b, d1, d2, and Area.
                            print(f"Found parallelogram: sides a={a}, b={b}; diagonals d1={d1}, d2={d2}; Area={area}")

    print(f"\nTotal number of distinct parallelograms found: {len(found_parallelograms)}")

solve()
<<<11>>>