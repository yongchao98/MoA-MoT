def solve_chesterton_query():
    """
    This function identifies and prints the name of the person
    who spoke a specific quote to G.K. Chesterton.
    """
    # The individual who made the statement to G.K. Chesterton
    speaker = "Lord Lloyd (George Lloyd, 1st Baron Lloyd)"

    # The location of the conversation
    location = "the Mount of Olives, Jerusalem"

    # The book where Chesterton recounts this event
    source = "G.K. Chesterton's Autobiography"

    # The quote itself
    quote = "'Well, anyhow, it must be obvious to anybody that the doctrine of the Fall is the only cheerful view of human life.'"

    # Printing the answer in a structured way
    print("According to his account in " + source + ", the person who spoke to G.K. Chesterton was:")
    print("Speaker: " + speaker)
    print("\nThey were standing together on " + location + " when Lord Lloyd said:")
    print("Quote: " + quote)


solve_chesterton_query()