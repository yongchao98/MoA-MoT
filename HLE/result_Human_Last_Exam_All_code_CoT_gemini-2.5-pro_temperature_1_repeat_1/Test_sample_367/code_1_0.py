import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def is_perfect_square(n):
    """Checks if a number is a perfect square and returns the root if it is."""
    if n < 0:
        return False, -1
    x = int(math.sqrt(n))
    return x * x == n, x

def count_parallelograms():
    """
    Counts the number of distinct parallelograms based on the given restrictions.
    A distinct parallelogram is defined by a unique set of {a, b, d1, d2}.
    """
    count = 0
    
    # Condition: 2a < a + b < 100, which implies a < b and a + b < 100.
    # As reasoned in the plan, a cannot be 1, so we start the loop from a = 2.
    # The maximum value for 'a' is 49, because if a=49, b can be 50 (since 49+50=99<100).
    for a in range(2, 50):
        # For a given 'a', b must be greater than 'a' and a + b < 100.
        # So, b ranges from a + 1 up to 99 - a.
        for b in range(a + 1, 100 - a):
            
            # Condition: a and b are coprime.
            if gcd(a, b) != 1:
                continue

            # From the parallelogram law: d1^2 + d2^2 = 2 * (a^2 + b^2)
            N = 2 * (a**2 + b**2)
            
            # Search for integer solutions (d1, d2). We assume d1 <= d2 to avoid duplicates.
            # This implies d1^2 <= N / 2, so d1 <= sqrt(N / 2).
            d1_limit = int(math.sqrt(N / 2))
            
            for d1 in range(1, d1_limit + 1):
                d2_squared = N - d1**2
                is_sq, d2 = is_perfect_square(d2_squared)
                
                if is_sq:
                    # Found a pair of integer diagonals (d1, d2). Now check remaining conditions.
                    
                    # Condition: Not a rectangle (diagonals are not equal).
                    if d1 == d2:
                        continue
                        
                    # Condition: The triangles (a, b, d1) and (a, b, d2) must be non-degenerate.
                    # This requires the strict triangle inequality: |a-b| < diagonal < a+b.
                    lower_bound = b - a
                    upper_bound = a + b
                    
                    if (d1 > lower_bound and d1 < upper_bound and
                        d2 > lower_bound and d2 < upper_bound):
                        # This set {a, b, d1, d2} defines a valid, distinct parallelogram.
                        count += 1
                        
    return count

final_count = count_parallelograms()
print(final_count)
<<<46>>>