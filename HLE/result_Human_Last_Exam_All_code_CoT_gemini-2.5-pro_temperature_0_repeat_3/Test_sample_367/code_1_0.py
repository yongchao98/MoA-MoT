import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def is_square(n):
    """Checks if a non-negative number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def find_sum_of_two_squares(N):
    """Finds integer pairs (x, y) such that x^2 + y^2 = N, with x > y > 0."""
    pairs = []
    # We only need to check y up to sqrt(N/2)
    limit = int(math.sqrt(N / 2))
    for y in range(1, limit + 1):
        x_sq = N - y * y
        if is_square(x_sq):
            x = int(math.sqrt(x_sq))
            # Ensure x > y to get unique pairs of squares
            if x > y:
                pairs.append((x, y))
    return pairs

def solve_parallelogram_problem():
    """
    Finds the number of distinct parallelograms satisfying the given restrictions.
    """
    found_pairs = set()

    # 1. Iterate through all possible side lengths a and b.
    # The condition 2a < a + b < 100 implies a < b and a + b < 100.
    for a in range(1, 50):
        for b in range(a + 1, 100):
            if a + b >= 100:
                break
            
            # 2. Sides a and b must be coprime.
            if gcd(a, b) != 1:
                continue

            # 3. From the parallelogram law, d1^2 + d2^2 = 2(a^2 + b^2).
            N = 2 * (a*a + b*b)
            
            # 4. Find all possible integer diagonal pairs (d1, d2).
            diag_pairs = find_sum_of_two_squares(N)
            
            for d1, d2 in diag_pairs:
                # A rectangle has equal diagonals. This is excluded.
                if d1 == d2:
                    continue

                # 5. Calculate K = ab*cos(theta) = (d1^2 - d2^2) / 4.
                # d1 and d2 must have the same parity for K to be an integer.
                # This is guaranteed by the nature of sums of two squares.
                K_numerator = d1*d1 - d2*d2
                if K_numerator % 4 != 0:
                    continue
                K = K_numerator // 4
                
                P = a * b
                
                # 6. Check for non-degenerate parallelogram: |cos(theta)| < 1 => |K| < P.
                if abs(K) >= P:
                    continue
                    
                # 7. Check for integer area: Area^2 = (ab)^2 - K^2 must be a perfect square.
                area_sq = P*P - K*K
                if is_square(area_sq):
                    # Found a valid parallelogram. Add the pair (a,b) to our set.
                    pair = tuple(sorted((a, b)))
                    found_pairs.add(pair)

    print(len(found_pairs))

solve_parallelogram_problem()