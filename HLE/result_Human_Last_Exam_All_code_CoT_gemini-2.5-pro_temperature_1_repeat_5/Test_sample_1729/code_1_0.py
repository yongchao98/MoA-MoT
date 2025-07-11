import sympy

def find_pm():
    """
    This function calculates the probability P_m symbolically.
    P_m is the probability that after removing two terms a_i and a_j from an
    arithmetic sequence of length 4m+2, the remaining 4m terms can be
    partitioned into m arithmetic sequences of 4 terms each.
    """

    # Define m as a symbolic variable
    m = sympy.Symbol('m', integer=True, positive=True)

    print("The probability P_m is calculated as the ratio of successful pairs (i,j) to the total possible pairs (i,j).")
    print("\nStep 1: Counting the number of successful pairs (i,j).")

    # A pair (i,j) is successful if the remaining 4m indices can be partitioned into m arithmetic progressions of length 4.
    # It can be shown that this is possible if the problem is viewed as partitioning the set of indices {1, ..., 4m+2}.
    # This leads to a combinatorial problem: arranging m blocks of 4 consecutive integers and 2 blocks of 1 integer (the removed i and j)
    # to form the full set of indices {1, 2, ..., 4m+2}.
    # The number of ways to do this is to choose 2 positions for the single-integer blocks out of m+2 total block positions.
    
    num_successful_pairs_expr = sympy.binomial(m + 2, 2)
    num_successful_pairs_factored = sympy.factor(num_successful_pairs_expr)
    
    print("The number of successful pairs is given by the binomial coefficient C(m+2, 2).")
    print(f"Number of successful pairs = {num_successful_pairs_expr} = {num_successful_pairs_factored}")
    
    print("\nStep 2: Counting the total number of pairs (i,j).")
    
    # The total number of pairs (i,j) that can be chosen from 4m+2 indices is given by the binomial coefficient C(4m+2, 2).
    total_pairs_expr = sympy.binomial(4*m + 2, 2)
    total_pairs_factored = sympy.factor(total_pairs_expr)
    
    print(f"The total number of pairs is C(4m+2, 2).")
    print(f"Total pairs = {total_pairs_expr} = {total_pairs_factored}")
    
    print("\nStep 3: Calculating the probability P_m.")
    
    # The probability P_m is the ratio of the two quantities.
    Pm = num_successful_pairs_expr / total_pairs_expr
    simplified_Pm = sympy.simplify(Pm)
    
    print("P_m = (Number of successful pairs) / (Total pairs)")
    print(f"P_m = ({num_successful_pairs_factored}) / ({total_pairs_factored})")

    # Display the simplified fraction and its components.
    
    numerator, denominator = simplified_Pm.as_numer_denom()
    
    print("\nThe simplified final formula for P_m is:")
    print(f"P_m = {simplified_Pm}")
    
    print("\nBreaking down the final equation:")
    print(f"Numerator: {numerator}")
    num_factors = sympy.factor_list(numerator)[1]
    # To handle the case where numerator is 1
    if not num_factors:
        print("Numerator factors: 1")
    else:
        print("Numerator factors:", ' * '.join([f"({factor**exp})" for factor, exp in num_factors]))
    
    print(f"Denominator: {denominator}")
    den_factors = sympy.factor_list(denominator)[1]
    print("Denominator factors:", ' * '.join([f"({factor**exp})" for factor, exp in den_factors]))


if __name__ == '__main__':
    find_pm()