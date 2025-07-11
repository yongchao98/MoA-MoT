def solve():
    """
    This function determines the query complexity for sorting bitstrings in two different regimes.
    The complexity is expressed in the (a,b,c) notation representing Theta(sqrt(N^a * (log N)^b * (log log N)^c)).

    Our analysis concluded that for both regimes, the optimal query complexity is Theta(N log N).
    Theta(N log N) = Theta(sqrt( (N log N)^2 )) = Theta(sqrt( N^2 * (log N)^2 * (log log N)^0 ))
    This corresponds to (a,b,c) = (2,2,0).
    """

    # For Regime 1: N = 2^sqrt(L)
    # The optimal strategy (either comparison sort or radix sort) results in a complexity of Theta(N log N).
    # This translates to (a,b,c) = (2,2,0).
    complexity1 = (2, 2, 0)

    # For Regime 2: N = 2^((log2 L)^2)
    # The radix sort strategy is less efficient than a simple comparison sort.
    # The optimal complexity is therefore that of comparison sort, which is Theta(N log N).
    # This translates to (a,b,c) = (2,2,0).
    complexity2 = (2, 2, 0)

    # The final answer is the combination of the two, in the specified format.
    # The problem asks for the answer in the format "(a,b,c),(a',b',c')"
    final_answer_string = f"({complexity1[0]},{complexity1[1]},{complexity1[2]}),({complexity2[0]},{complexity2[1]},{complexity2[2]})"
    
    print(final_answer_string)
    print("<<<" + final_answer_string + ">>>")

solve()