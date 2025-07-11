import math

def get_prime_factors(n):
    """Returns the set of prime factors of an integer n."""
    factors = set()
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            factors.add(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.add(temp)
    return factors

def check_a_conditions(a):
    """
    Checks if 'a' is odd and all its prime factors are congruent to 1 mod 4.
    """
    if a % 2 == 0:
        return False
    prime_factors = get_prime_factors(a)
    for factor in prime_factors:
        if factor % 4 != 1:
            return False
    return True

def get_crossing_number(p, q):
    """
    Calculates the crossing number of the 2-bridge knot K(p,q)
    by summing the coefficients of the continued fraction expansion of p/q.
    """
    if q == 0:
        return 0
    # Normalize q to be the smaller of q and p-q for calculation, as it doesn't change the knot type.
    # Note: K(p,q) and K(p, p-q) can have different CF expansions but must have the same crossing number.
    # Standard practice takes 0 < q < p/2.
    q_norm = q % p
    if q_norm > p // 2:
        q_norm = p - q_norm
        
    temp_p, temp_q = p, q_norm
    c_num = 0
    while temp_q > 0:
        c_num += temp_p // temp_q
        temp_p, temp_q = temp_q, temp_p % temp_q
    return c_num

def find_knots():
    """
    Finds the 2-bridge knots with crossing number <= 13 satisfying the conditions.
    """
    max_crossing_number = 13
    found_knots_params = []
    
    # Start with the unknot (c=0)
    # The theorem is for non-trivial knots, so we add the unknot manually.
    print("Found knot: Unknot (0_1), crossing number 0")
    
    # We estimate an upper bound for 'a'. If c <= 13, p must be < ~8000, so a < 90.
    for a in range(3, 91):
        if not check_a_conditions(a):
            continue
        
        p = a * a
        # We need to find an even q such that q^2 = -1 (mod p) and gcd(p,q)=1
        # We can iterate through q and check
        # Since p=a^2 is odd, if we find a q, p-q will be odd. So we don't double count.
        for q in range(2, p, 2): # q must be even
            if math.gcd(p, q) != 1:
                continue

            if (q * q) % p == (p - 1):
                # Found a candidate knot K(p,q)
                crossing_number = get_crossing_number(p, q)
                
                if crossing_number <= max_crossing_number:
                    found_knots_params.append({'p': p, 'q': q, 'c': crossing_number})
                    print(f"Found knot: K({p}, {q}), crossing number {crossing_number}")
                # As 'a' grows, 'p' grows, and 'c' tends to grow. We can break early
                # but searching the full range up to a=90 is fast enough and safer.
                # Since q solutions are sparse, we don't need to optimize this inner loop heavily.
    
    # The list now contains parameters of the knots.
    # From our analysis, each 'a' value gives a unique knot, so the list doesn't contain duplicates.
    count = 1 + len(found_knots_params)
    print("\nTotal number of such knots found.")
    
    # The question asked to output the numbers in the final equation.
    equation_parts = ["1 (unknot)"]
    for _ in found_knots_params:
        equation_parts.append("1")
    
    print(f"{' + '.join(equation_parts)} = {count}")
    return count

if __name__ == '__main__':
    total_knots = find_knots()
    # The final answer format is specific, but not requested for the main script block.
    # The final print within find_knots serves the purpose.
    # print(f"\nFinal Answer: {total_knots}")
