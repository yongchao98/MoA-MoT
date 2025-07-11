import math

def count_knots(max_crossing_number):
    """
    Calculates the number of 2-bridge knots with crossing number up to max_crossing_number
    that admit two disjoint non-parallel embedded minimal genus Seifert surfaces.
    """
    total_knots = 0
    contributions = []

    # The crossing number 'c' must be even.
    # We iterate through even crossing numbers from 2 up to max_crossing_number.
    for c in range(2, max_crossing_number + 1, 2):
        knots_for_c = 0
        k = c // 2

        # The number of terms in the continued fraction 'n' must be even.
        # Also, since each part b_i >= 1, the sum k must be at least n.
        for n in range(2, k + 1, 2):
            # N_compositions is the number of ways to write k as a sum of n positive integers.
            # This is given by the combination formula C(k-1, n-1).
            try:
                N_compositions = math.comb(k - 1, n - 1)
            except ValueError:
                N_compositions = 0
            
            # N_palindromes is the number of those compositions that are palindromic.
            # A composition (b_1, ..., b_n) is palindromic if b_i = b_{n-i+1}.
            # For n even, this requires 2 * sum(b_1 to b_{n/2}) = k.
            # This is only possible if k is even.
            N_palindromes = 0
            if k % 2 == 0:
                # The number of palindromic compositions is the number of compositions of k/2 into n/2 parts.
                try:
                    N_palindromes = math.comb(k // 2 - 1, n // 2 - 1)
                except ValueError:
                    N_palindromes = 0
            
            # The number of distinct knots is the number of compositions up to reversal.
            # This is (number of non-palindromic + 2*number of palindromic)/2
            # = (total compositions + palindromic compositions) / 2
            knots_for_cn = (N_compositions + N_palindromes) // 2
            knots_for_c += knots_for_cn
            
        print(f"For crossing number {c}, there are {knots_for_c} such knots.")
        total_knots += knots_for_c
        contributions.append(str(knots_for_c))

    print("\nThe final calculation is:")
    final_equation = " + ".join(contributions) + f" = {total_knots}"
    print(final_equation)
    print(f"\nThe total number of such knots with crossing number at most {max_crossing_number} is: {total_knots}")
    
    return total_knots

if __name__ == '__main__':
    MAX_CROSSING = 13
    final_answer = count_knots(MAX_CROSSING)
    print(f"<<<{final_answer}>>>")
