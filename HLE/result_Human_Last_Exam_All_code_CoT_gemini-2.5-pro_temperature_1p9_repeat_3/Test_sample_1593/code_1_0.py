def solve():
    """
    This function prints the determined query complexities in (a,b,c) notation.
    """
    # For the regime N = 2^sqrt(L), the complexity is Theta(N log N).
    # In abc-notation, this is (a=2, b=2, c=0).
    regime_1_abc = (2, 2, 0)

    # For the regime N = 2^((log_2 L)^2), the complexity is also Theta(N log N).
    # In abc-notation, this is (a=2, b=2, c=0).
    regime_2_abc = (2, 2, 0)
    
    # Format the output string as specified in the example "(2,0,0),(2,1,-1)"
    # The prompt also says to "output each number in the final equation",
    # which we interpret as printing the components of the (a,b,c) tuples.
    answer_string = f"({regime_1_abc[0]},{regime_1_abc[1]},{regime_1_abc[2]}),({regime_2_abc[0]},{regime_2_abc[1]},{regime_2_abc[2]})"
    
    print(answer_string)

solve()