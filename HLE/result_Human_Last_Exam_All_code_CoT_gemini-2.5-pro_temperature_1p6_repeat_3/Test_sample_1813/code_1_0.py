import math
from fractions import Fraction

def get_continued_fraction_coeffs(r: Fraction):
    """
    Computes the list of coefficients for the continued fraction expansion of a rational number.
    Example: 4/7 -> [0, 1, 1, 3]
    """
    coeffs = []
    num = r.numerator
    den = r.denominator
    
    while den != 0:
        quotient = num // den
        remainder = num % den
        coeffs.append(quotient)
        num = den
        den = remainder
        
    return coeffs

def main():
    """
    Main function to compute and display the continued fraction associated
    with the generalized Markov number m_{4/7}.
    """
    # The rational number in question is 4/7.
    r = Fraction(4, 7)
    
    # Step 1: Compute the continued fraction of r = 4/7.
    # The coefficients will be of the form [a_0; a_1, a_2, ..., a_n].
    all_coeffs = get_continued_fraction_coeffs(r)
    
    # Since 0 < 4/7 < 1, the first coefficient a_0 is 0.
    # We are interested in the block of coefficients A = (a_1, ..., a_n).
    A = all_coeffs[1:]
    
    # Step 2: Construct the reverse of block A.
    A_rev = A[::-1]
    
    # Step 3: The repeating block of the associated continued fraction is A followed by A_rev.
    repeating_block = A + A_rev
    
    # Step 4: Format and print the final result.
    # The associated continued fraction is [0; overline{repeating_block}].
    # We will print this in a readable format, showing all the numbers.
    
    print(f"The continued fraction expansion of {r} is: [{all_coeffs[0]}; {', '.join(map(str, A))}]")
    
    repeating_block_str = ', '.join(map(str, repeating_block))
    
    print("\nThe continued fraction associated with the generalized Markov number m_{4/7} is:")
    # We represent the repeating nature by indicating the block.
    # The final equation form is [0; a, b, c, ...], where the block [a, b, c, ...] repeats.
    print(f"m_{{4/7}} = [0; {repeating_block_str}, ...]")
    print(f"This is a periodic continued fraction where the block [{repeating_block_str}] repeats infinitely.")
    print(f"In standard notation: m_{{4/7}} = [0; \\overline{{{repeating_block_str}}}]")


if __name__ == "__main__":
    main()
