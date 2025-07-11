import decimal

def solve():
    """
    Calculates the number of sequences and the parameters (a1, b1) for the largest term
    to represent 10^100 as a sum of tetration-power terms.
    """
    # Set a high precision for decimal calculations to handle logarithms of large numbers accurately.
    decimal.getcontext().prec = 150

    LOG2_CONST = decimal.Decimal(2).ln()

    def log2_dec(x):
        """Calculates the base-2 logarithm of a decimal number."""
        if x <= 0:
            return decimal.Decimal('-Infinity')
        return x.ln() / LOG2_CONST

    # Start with the log2 of the number N = 10^100.
    # log2(10^100) = 100 * log2(10)
    log2_N = decimal.Decimal(100) * log2_dec(decimal.Decimal(10))

    count = 0
    a1, b1 = 0, 0

    # Loop until the remaining number is negligible (i.e., less than 1, so its log2 < 0).
    while log2_N > 0:
        count += 1

        # Find the largest 'a' such that tet(2, a) <= N_rem.
        # This is equivalent to finding the largest integer 'a' such that:
        # a <= 1 + log2(log2(N_rem))
        # We need log2_N to be > 1 for log2(log2_N) to be positive. If log2_N < 1,
        # no tet(2,a) with a>=1 can be a part of the term.
        if log2_N <= 1:
            break
        
        log2_log2_N = log2_dec(log2_N)
        if log2_log2_N.is_infinite(): # This happens if log2_N is too small.
            break
            
        a = int(1 + log2_log2_N)

        if a < 1:
            break
            
        # The log2 of the tetration part of the term: log2(tet(2, a)) = 2^(a-1)
        log2_tet_a = decimal.Decimal(2)**(a - 1)
        
        # Find the largest integer 'b' such that tet(2, a) * 2^b <= N_rem.
        # This implies b <= log2_N - log2(tet(2, a)).
        # We assume b must be non-negative.
        if log2_N < log2_tet_a:
             # This should not be reached due to how 'a' is chosen.
             continue
        b = int(log2_N - log2_tet_a)

        # The constraint is tet(2, a) > 2^b, which means log2(tet(2, a)) > b.
        # The greedy choice of a and b ensures this.

        # Store the a and b for the first (and largest) term.
        if count == 1:
            a1 = a
            b1 = b

        # Update log2_N for the remainder using a stable logarithmic formula.
        # N_new = N_old - term = (term * 2^f) - term = term * (2^f - 1)
        # where f is the fractional part of the logarithm.
        f = log2_N - log2_tet_a - decimal.Decimal(b)
        
        # The factor by which the term was overestimated.
        term_in_log = (decimal.Decimal(2)**f) - 1
        
        if term_in_log <= 0:
            # Remainder is zero or negative (due to precision), so we are done.
            break
            
        log2_rem_factor = log2_dec(term_in_log)
        
        # The log of the new remainder.
        log2_N = log2_tet_a + decimal.Decimal(b) + log2_rem_factor

    print(f"{count} {a1} {b1}")

solve()
<<<110 9 76>>>