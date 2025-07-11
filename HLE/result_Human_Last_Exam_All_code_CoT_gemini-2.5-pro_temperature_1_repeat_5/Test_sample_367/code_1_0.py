import math

def main():
    """
    This script calculates the number of distinct parallelograms based on the given restrictions.
    """

    def gcd(a, b):
        """Computes the greatest common divisor of a and b."""
        while b:
            a, b = b, a % b
        return a

    def is_perfect_square(n):
        """Checks if n is a perfect square and returns (bool, sqrt(n))."""
        if n < 0:
            return False, -1
        if n == 0:
            return True, 0
        x = int(math.sqrt(n))
        return x * x == n, x

    parallelogram_count = 0

    # Condition 2: Iterate through coprime integers a, b with 2a < a + b < 100.
    # This is equivalent to 1 <= a < b and a + b < 100.
    for b in range(2, 100):
        for a in range(1, b):
            # Check the sum constraint
            if a + b >= 100:
                break
            
            # Check the coprime constraint
            if gcd(a, b) != 1:
                continue

            # For each valid side pair (a, b), find possible integer diagonals.
            # Condition 4: d1 and d2 are integers.
            # From the parallelogram law: d1^2 + d2^2 = 2*(a^2 + b^2).
            N = 2 * (a*a + b*b)
            
            # We search for integer solutions {d1, d2} for the equation x^2 + y^2 = N.
            # We assume d1 <= d2 to count unique pairs of diagonals.
            limit = int(math.sqrt(N / 2))
            for d1 in range(1, limit + 1):
                d2_squared = N - d1*d1
                is_sq, d2 = is_perfect_square(d2_squared)
                
                if is_sq:
                    # We have found an integer diagonal pair {d1, d2}.
                    
                    # A valid parallelogram can only be formed if d2 <= a+b.
                    # This single check implies |a-b| <= d1 as well.
                    if d2 > a + b:
                        continue
                    
                    # Condition 1: The parallelogram is not a rectangle.
                    # A rectangle has equal diagonals (d1 == d2).
                    if d1 == d2:
                        continue
                    
                    # If all conditions are met, we have found a distinct parallelogram.
                    parallelogram_count += 1
    
    # The final equation is the total count.
    print(f"{parallelogram_count}")

if __name__ == "__main__":
    main()