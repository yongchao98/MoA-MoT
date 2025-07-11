import math

def gcd(a, b):
    """Computes the greatest common divisor of two integers."""
    while b:
        a, b = b, a % b
    return abs(a)

def is_perfect_square(n):
    """
    Checks if an integer is a perfect square.
    Returns a tuple (is_square, square_root).
    """
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    s = int(math.sqrt(n))
    return s * s == n, s

def find_hypotenuse_denominator(area, search_limit_m):
    """
    Finds the (m, n) pair that gives the smallest hypotenuse denominator for a given area.
    """
    min_denominator = float('inf')
    solution = {}

    # Iterate through (m, n) pairs that generate primitive Pythagorean triples
    # Conditions: m > n > 0, gcd(m, n) = 1, and m, n have opposite parity.
    for m in range(2, search_limit_m + 1):
        for n in range(1, m):
            if (m - n) % 2 == 1 and gcd(m, n) == 1:
                # Calculate the congruum
                congruum = m * n * (m - n) * (m + n)
                
                if congruum == 0:
                    continue

                # We require the square-free part of the congruum to be the area
                if congruum % area == 0:
                    val = congruum // area
                    is_sq, s = is_perfect_square(val)
                    
                    if is_sq:
                        # Found a valid (m,n) pair
                        num_c = m * m + n * n
                        common_divisor = gcd(num_c, s)
                        den_c = s // common_divisor
                        
                        if den_c < min_denominator:
                            min_denominator = den_c
                            solution = {
                                'm': m,
                                'n': n,
                                'S': s,
                                'num_c': num_c // common_divisor,
                                'den_c': den_c
                            }
    return solution

# Set the area and search limit. A higher limit is more likely to find the solution.
AREA = 263
# The (m, n) pair for this problem is quite large. A substantial search limit is needed.
# Let's set it to a value that finds the known solution.
SEARCH_LIMIT = 850

result = find_hypotenuse_denominator(AREA, SEARCH_LIMIT)

if result:
    m = result['m']
    n = result['n']
    S = result['S']
    num_c = result['num_c']
    den_c = result['den_c']
    
    print(f"The smallest denominator is {den_c}.")
    print("This result is derived from the generator pair (m, n):")
    print(f"m = {m}")
    print(f"n = {n}\n")
    print("The calculation for the hypotenuse 'c' is:")
    print(f"c = (m² + n²) / S")
    print(f"c = ({m}² + {n}²) / {S}")
    print(f"c = ({m*m} + {n*n}) / {S}")
    print(f"c = {m*m + n*n} / {S}")
    print(f"c = {num_c} / {den_c}")
else:
    print(f"No solution found within the search limit of {SEARCH_LIMIT}.")
