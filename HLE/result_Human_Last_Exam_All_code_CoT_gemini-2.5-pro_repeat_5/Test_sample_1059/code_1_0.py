def solve():
    """
    Determines which of the given mathematical sets are necessarily closepact.

    The analysis of each option is as follows:
    A. The set of real numbers: No.
    B. The set of integers: No.
    C. A finite subset of the complex numbers: Yes.
    D. The set of all 1/n where n is a nonzero integer: No.
    E. The set containing a Cauchy sequence in the rationals: No.
    F. The set containing a bounded monotonic sequence in the real numbers: No.
    G. The set containing a bounded monotonic sequence and its limit point in the real numbers: No.
    H. The set containing a positive real sequence and its limit point: No.
    I. An open interval in the reals: No.
    J. A closed interval in the reals: Yes.
    K. A bounded measurable subset of the real numbers: No.
    L. A bounded non-measurable subset of the real numbers: No.
    M. The Cantor Set: No.

    The correct choices are C and J.
    """
    answer = "CJ"
    print(answer)

solve()