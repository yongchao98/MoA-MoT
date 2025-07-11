import math

def find_parallelograms():
    """
    Finds and counts the number of distinct parallelograms based on the problem's restrictions.
    """
    count = 0
    print("Found parallelograms (a, b, d1, d2):")
    # Loop through possible side lengths a and b.
    # a < b and a + b < 100.
    # b can go up to 98 (if b=98, a=1, a+b=99).
    for b in range(2, 99):
        # a must be less than b and a < 100-b.
        limit_a = min(b, 100 - b)
        for a in range(1, limit_a):
            # Condition 2: a and b are coprime.
            if math.gcd(a, b) != 1:
                continue

            # Condition 1: Parallelogram is not a rectangle.
            # This means a^2 + b^2 is not a perfect square.
            c_sq = a**2 + b**2
            c = int(math.sqrt(c_sq))
            if c**2 == c_sq:
                continue

            # From the parallelogram law, 2*(a^2 + b^2) = d1^2 + d2^2.
            N = 2 * c_sq
            
            # Find all pairs of integer diagonals (d1, d2) for the given a, b.
            # We iterate d1 up to sqrt(N/2) to find pairs (d1, d2) where d1 <= d2.
            limit_d1 = int(math.sqrt(N // 2))
            for d1 in range(1, limit_d1 + 1):
                d2_sq = N - d1**2
                d2 = int(math.sqrt(d2_sq))

                if d2**2 == d2_sq:
                    # We found a pair of integer diagonals (d1, d2).
                    
                    # Condition: Area must be a non-zero integer.
                    # This implies two sub-conditions:
                    # 1. The parallelogram is not degenerate (strict triangle inequality).
                    #    |a-b| < d1 < a+b and |a-b| < d2 < a+b
                    b_minus_a = b - a
                    a_plus_b = a + b
                    if not (d1 > b_minus_a and d2 < a_plus_b):
                        continue

                    # 2. The area is an integer, not a half-integer. This holds if
                    #    d1 and (a+b) have the same parity. d2 will automatically have
                    #    the same parity as d1 if 2(a^2+b^2) = d1^2+d2^2.
                    if (d1 % 2) != (a_plus_b % 2):
                        continue
                    
                    # All conditions are met. This is a valid parallelogram.
                    # Print the numbers in the final equation: d1^2 + d2^2 = 2*(a^2+b^2)
                    print(f"a={a}, b={b}, d1={d1}, d2={d2}. Equation: {d1}^2 + {d2}^2 = 2*({a}^2 + {b}^2) => {d1**2 + d2**2} = {N}")
                    count += 1
    
    print("\nTotal number of distinct parallelograms:")
    print(count)

if __name__ == '__main__':
    find_parallelograms()
<<<112>>>