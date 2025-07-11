import math

def solve():
    """
    Calculates the sum over all natural dimensions of the specified ratio.
    """
    # The constant K is the expected distance E[sqrt(c_i^2 + c_j^2)] for c_i, c_j ~ U(-1,1).
    # Its value is K = (1/3) * (sqrt(2) + log(1 + sqrt(2))).
    # math.asinh(1) is a precise way to compute log(1 + sqrt(2)).
    K = (1.0/3.0) * (math.sqrt(2) + math.asinh(1))

    total_sum = 0.0
    d = 1
    
    # We sum the terms R_d of the series until they become too small to affect the result,
    # ensuring the sum is accurate to the required precision.
    while True:
        # The formula for the term R_d in the series is:
        # R_d = (d+1) * (0.5)^d / (d! * (1 + (d-1)*K))
        try:
            term_numerator = (d + 1) * (0.5)**d
            term_denominator = math.factorial(d) * (1 + (d - 1) * K)
            term = term_numerator / term_denominator
        except OverflowError:
            # This is a safeguard; the loop will break due to the 'term' size
            # long before math.factorial overflows.
            break
            
        # When a term becomes sufficiently small, its contribution to the sum is negligible.
        if term < 1e-12:
            break
            
        total_sum += term
        d += 1

    # The problem asks for the final sum determined with three-decimal precision.
    # The numbers in the final equation are the terms of the series (R_1, R_2, ...)
    # and their sum (S). The code calculates this sum S.
    print(f"{total_sum:.3f}")

solve()