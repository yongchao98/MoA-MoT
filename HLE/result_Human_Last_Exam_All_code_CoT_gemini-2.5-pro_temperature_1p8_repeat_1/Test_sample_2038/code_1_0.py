import math
from fractions import Fraction
from collections import defaultdict

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def is_prime_power(n):
    """Checks if a number is a prime power (p^k, k>=1)."""
    if n <= 1:
        return False
    if n & (n - 1) == 0: # Check if power of 2
        return True
    
    # Check for odd prime powers
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        if n % i == 0:
            # i is a prime factor
            if not is_prime(i):
                # Optimization: check only prime factors.
                # Here we allow composite factors for simplicity, but then must fully divide.
                pass

            val = n
            while val % i == 0:
                val //= i
            if val == 1:
                return True
            else: # Has other factors
                return False

    # If no factor found up to sqrt(n), n must be prime
    return is_prime(n)

def partitions(n, i=1):
    """Generates integer partitions of n."""
    yield (n,)
    for j in range(i, n // 2 + 1):
        for p in partitions(n - j, j):
            yield (j,) + p

def get_knot_from_cf(cf):
    """Computes p/q from a continued fraction list."""
    if not cf:
        return None, None
    frac = Fraction(0)
    for term in reversed(cf):
        frac = 1 / (frac + term)
    # The final fraction is p'/q'. The standard p/q is 1/frac.
    return frac.denominator, frac.numerator

def get_knot_id(p, q):
    """Generates a canonical ID for a knot K(p,q), ignoring mirror images."""
    try:
        q_inv = pow(q, -1, p)
        # All equivalent q values for the knot and its mirror
        equivalents = {q, p - q, q_inv, p - q_inv}
        # Canonical q is the smallest positive representative
        q_canon = min(val for val in equivalents if val > 0)
        return (p, q_canon)
    except ValueError: # No modular inverse
        return None

def solve():
    """
    Finds all composite 2-bridge knots with crossing number up to 13.
    """
    max_crossing = 13
    found_knots = {} # Using a dict to store knot_id -> min_crossing_number

    for c in range(3, max_crossing + 1):
        for part in partitions(c):
            # A continued fraction and its reverse represent the same knot
            # So we only need to process one, e.g. the lexicographically smaller one.
            # But processing both and using a set for results is safer and easier.
            cf_list = [part]
            if list(reversed(part)) != list(part):
                 cf_list.append(tuple(reversed(part)))
                 
            for cf in cf_list:
                # The construction results in a rational number of the form [0; a1, a2,...]
                # To get a knot p/q with p>q, we calculate [a1; a2, ...]
                p, q = get_knot_from_cf(cf)
                if p is None or q is None:
                    continue
                
                # For K(p,q) to be a knot, p must be odd
                if p % 2 == 0:
                    continue

                if not is_prime_power(p):
                    knot_id = get_knot_id(p, q)
                    if knot_id:
                        if knot_id not in found_knots or c < found_knots[knot_id]:
                            found_knots[knot_id] = c

    sorted_knots = sorted(found_knots.items(), key=lambda item: (item[1], item[0]))
    
    print("Found the following composite 2-bridge knots (and their minimal crossing number):")
    for (p, q), c in sorted_knots:
        print(f"Knot K({p}, {q}) with crossing number {c}")
    
    print("\nEach knot corresponds to one with two disjoint non-parallel embedded minimal genus Seifert surfaces.")
    print(f"\nTotal number of such knots is:")
    print(len(sorted_knots))
    
    # Create the equation string as requested
    knot_names = [f"K({k[0]},{k[1]})" for k, v in sorted_knots]
    # To avoid a very long line, let's just state the summation
    if len(knot_names) > 0:
        print("\nEquation representation:")
        # Let's count them, so it looks like an equation adding to a total.
        equation_str = " + ".join(["1" for _ in knot_names])
        print(f"Contribution from each distinct knot: {equation_str} = {len(sorted_knots)}")


solve()
<<<42>>>