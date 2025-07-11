import sys

def solve():
    """
    This function prints the derived formula for P(n).
    
    The derivation steps are outlined in the thought process. It involves:
    1. Taking the logarithm of Q(n).
    2. Using the Euler-Maclaurin formula to find the asymptotic expansion of ln(Q(n)).
    3. Defining delta(n) as the variable part of the expansion.
    4. Expanding exp(delta(n)) as a series.
    5. Identifying the terms that must be included in P(n) to achieve the desired error bound.
    
    The terms of order n^-2 and n^-3 in the expansion of exp(delta(n)) beyond the already included (1 + L/(2n)) are:
    P_2 = (3*L**2 - 2*L + 2) / (24*n**2)
    P_3 = (L**3 - 2*L**2 + 2*L) / (48*n**3)
    
    P(n) is the sum of these two terms.
    """
    
    # Define the numerator and denominator for the n^-2 term of P(n)
    term1_numerator = "3*L**2 - 2*L + 2"
    term1_denominator = "24*n**2"
    
    # Define the numerator and denominator for the n^-3 term of P(n)
    term2_numerator = "L**3 - 2*L**2 + 2*L"
    term2_denominator = "48*n**3"

    # Construct and print the final formula for P(n)
    # The numbers in the equation are explicitly written out as requested.
    # P(n) is the sum of the term of order n^{-2} and the term of order n^{-3}.
    # The term of order n^{-2} is (3*L^2 - 2*L + 2)/(24*n^2)
    # The term of order n^{-3} is (L^3 - 2*L^2 + 2*L)/(48*n^3)
    
    final_formula = f"P(n) = ({term1_numerator}) / ({term1_denominator}) + ({term2_numerator}) / ({term2_denominator})"
    
    print(final_formula)

solve()