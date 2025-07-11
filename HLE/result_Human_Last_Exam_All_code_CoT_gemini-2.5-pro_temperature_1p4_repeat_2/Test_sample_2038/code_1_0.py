import math

def get_continued_fraction_coeffs(n, d):
    """Computes the coefficients of the continued fraction for n/d."""
    coeffs = []
    while d != 0:
        a = n // d
        coeffs.append(a)
        n, d = d, n % d
    return coeffs

def sum_cf_coeffs(p, q):
    """Calculates the sum of the continued fraction coefficients."""
    coeffs = get_continued_fraction_coeffs(p, q)
    return sum(coeffs)

def modular_inverse(a, m):
    """Computes the modular inverse of a modulo m using Python's pow function."""
    try:
        # Requires Python 3.8+
        return pow(a, -1, m)
    except ValueError:
        # This will be raised if the inverse does not exist.
        return None

def get_canonical_representative(p, q):
    """
    Finds the canonical representative for the knot K(p,q), treating mirrors as identical.
    The equivalence class includes {q, q_inv, p-q, (p-q)_inv}.
    """
    q_inv = modular_inverse(q, p)
    if q_inv is None: return None

    q_mirror = p - q
    q_mirror_inv = modular_inverse(q_mirror, p)
    if q_mirror_inv is None: return None

    return (p, min(q, q_inv, q_mirror, q_mirror_inv))

def solve():
    """
    Finds the number of 2-bridge knots with crossing number (approximated by CF sum)
    at most 13 that are fibered.
    """
    max_crossing_sum = 13
    fibered_knots = set()

    # The largest p for a crossing number sum of 13 corresponds to the
    # continued fraction [1]*13, giving p=377 (a Fibonacci number).
    max_p = 377

    # Iterate through all possible p (odd integers) and q.
    for p in range(3, max_p + 1, 2):
        for q in range(1, p):
            if math.gcd(p, q) != 1:
                continue

            # We use the sum of continued fraction coefficients as a proxy for the crossing number.
            cf_sum = sum_cf_coeffs(p, q)

            if cf_sum <= max_crossing_sum:
                # A 2-bridge knot K(p,q) is fibered if q^2 is congruent to +/- 1 (mod p).
                q_squared_mod_p = (q * q) % p
                if q_squared_mod_p == 1 or q_squared_mod_p == p - 1:
                    # Add the knot's canonical representative to a set to count uniques.
                    canonical_rep = get_canonical_representative(p, q)
                    if canonical_rep:
                        fibered_knots.add(canonical_rep)
    
    # Print the result
    knot_list = sorted(list(fibered_knots))
    print(f"Found {len(knot_list)} unique 2-bridge knots that satisfy the conditions.")
    print("The canonical representatives (p,q) for these knots are:")
    for knot in knot_list:
      print(knot)
    
    total = len(knot_list)
    equation = " + ".join(["1"] * total)
    print("\nIn equation form, the count is:")
    print(f"{equation} = {total}")
    
    return total

if __name__ == '__main__':
    final_answer = solve()
    # The final answer in the required format
    # print(f"<<<{final_answer}>>>")
    # This part is for the platform, the user will see the print outputs from solve()

solve()
<<<30>>>