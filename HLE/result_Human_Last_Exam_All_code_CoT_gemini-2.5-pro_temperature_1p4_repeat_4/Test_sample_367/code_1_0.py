import math

def find_sum_of_two_squares(n, parity=None):
    """
    Finds all pairs (x, y) with x > y > 0 such that x^2 + y^2 = n.
    'parity' can be 'odd' or 'even' to filter x and y.
    """
    reps = []
    limit = int(math.sqrt(n / 2))
    for y in range(1, limit + 1):
        if parity == 'odd' and y % 2 == 0:
            continue
        if parity == 'even' and y % 2 != 0:
            continue

        x_sq = n - y * y
        if x_sq <= 0:
            continue
        
        x = int(math.sqrt(x_sq))

        if x * x == x_sq:
            if parity == 'odd' and x % 2 == 0:
                continue
            if parity == 'even' and x % 2 != 0:
                continue
                
            if x > y:
                reps.append(tuple(sorted((x, y), reverse=True)))
    return list(set(reps))


def solve():
    """
    Finds the number of distinct parallelograms based on the given restrictions.
    """
    count = 0
    solutions = []

    # Case 1: a, b are both odd.
    for a in range(1, 50, 2):
        for b in range(a + 2, 100 - a, 2):
            if math.gcd(a, b) != 1:
                continue

            # Diagonals are even: d1=2k, d2=2l. a^2+b^2=2(k^2+l^2)
            X = (a * a + b * b) // 2
            sum_sq_reps = find_sum_of_two_squares(X)

            for k, l in sum_sq_reps:
                c1_sq = (k * k - l * l)
                ab = a * b
                
                # Check if (c1, A, ab) is a Pythagorean triple
                area_sq = ab * ab - c1_sq * c1_sq
                if area_sq > 0:
                    area = int(math.sqrt(area_sq))
                    if area * area == area_sq:
                        count += 1
                        solutions.append({'a': a, 'b': b, 'd1': 2 * k, 'd2': 2 * l})

    # Case 2: a is even, b is odd
    for a in range(2, 50, 2):
        for b in range(a + 1, 100 - a, 2):
            if math.gcd(a, b) != 1:
                continue
            
            # Diagonals d1, d2 are odd. d1^2+d2^2 = 2(a^2+b^2)
            Y = 2 * (a * a + b * b)
            sum_sq_reps_odd = find_sum_of_two_squares(Y, parity='odd')
            
            for d1, d2 in sum_sq_reps_odd:
                c1 = d1*d1 - d2*d2 # This is d1^2 - d2^2
                
                # (c1, 4A, 4ab) is a Pythagorean triple
                term1_sq = (4 * a * b)**2
                term2_sq = c1**2

                four_A_sq = term1_sq - term2_sq
                if four_A_sq > 0:
                    # Check if four_A_sq is a perfect square divisible by 16
                    if four_A_sq % 16 == 0:
                        A_sq = four_A_sq // 16
                        A = int(math.sqrt(A_sq))
                        if A * A == A_sq:
                           count += 1
                           solutions.append({'a': a, 'b': b, 'd1': d1, 'd2': d2})

    print("Found the following parallelograms (sides a, b and diagonals d1, d2):")
    for s in solutions:
        # Sort diagonals for consistent output
        diagonals = tuple(sorted((s['d1'], s['d2']), reverse=True))
        print(f"a={s['a']}, b={s['b']}; d1={diagonals[0]}, d2={diagonals[1]}")
    
    # The final "equation" requested is simply the sum of found parallelograms.
    equation_parts = [1] * count
    equation_str = " + ".join(map(str, equation_parts))
    print(f"\nThe final count is the result of the sum: {equation_str} = {count}")


solve()
<<<9>>>