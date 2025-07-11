def solve():
    """
    Based on the analysis that in a metric space X, a subset Y is 'closepact'
    if and only if it is compact, we evaluate each option. A subset of R, C, or Q
    is compact if and only if it is closed and bounded (or totally bounded for Q).
    """
    
    # A. The set of real numbers: Not bounded. Not compact.
    # B. The set of integers: Not bounded. Not compact.
    # C. A finite subset of the complex numbers: Always closed and bounded. Compact.
    # D. The set of all 1/n where n is a nonzero integer: Not closed (limit point 0 is missing). Not compact.
    # E. The set containing a Cauchy sequence in the rationals: The set {q_n} is compact in Q
    #    (it's totally bounded, and it's closed in Q). Compact.
    # F. The set containing a bounded monotonic sequence in R: Not necessarily closed (limit may be missing). Not compact.
    # G. The set containing a bounded monotonic sequence and its limit point in R: Closed and bounded. Compact.
    # H. The set containing a positive real sequence and its limit point: Closed and bounded. Compact.
    # I. An open interval in the reals: Not closed. Not compact.
    # J. A closed interval in the reals: Closed and bounded. Compact.
    # K. A bounded measurable subset of R: Not necessarily closed. Not compact.
    # L. A bounded non-measurable subset of R: Must not be closed. Not compact.
    # M. The Cantor Set: Closed and bounded. Compact.
    
    answer = "CEGHJM"
    print(answer)

solve()