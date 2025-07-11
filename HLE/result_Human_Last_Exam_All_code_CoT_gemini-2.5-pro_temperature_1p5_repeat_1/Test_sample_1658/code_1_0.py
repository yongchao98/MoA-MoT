def solve_riddle():
    """
    This script identifies the person described in the user's query
    by connecting the historical figures mentioned.
    """

    # Key figures from the historical context
    host_writer = "Léon Bloy"
    artist_in_question = "Georges Rouault"
    thinker = "Jacques Maritain"
    author_of_quote = "Raïssa Maritain" # Jacques Maritain's wife

    # Explanation of the connection
    print(f"The home of the writer {host_writer} was a gathering place for artists and intellectuals.")
    print(f"The artist described in the quote is {artist_in_question}.")
    print(f"The 'prominent figure of European thought' who visited Bloy's home between 1905 and 1909 and became a close friend of both men was the philosopher {thinker}.")
    print(f"The quote itself was written by his wife, {author_of_quote}, in her memoirs about their shared experiences.")
    print("\nTherefore, the thinker is:")
    print(thinker)

solve_riddle()