def find_quote_origin():
    """
    This function identifies and prints the person who spoke a specific
    quote to G.K. Chesterton.
    """
    # The person who spoke the words to G.K. Chesterton.
    speaker = "his wife, Frances Chesterton"

    # Contextual information
    recipient = "G.K. Chesterton"
    location = "the Mount of Olives"
    quote = "'Well, anyhow, it must be obvious to anybody that the doctrine of the Fall is the only cheerful view of human life.'"

    # The final answer statement
    answer_statement = (
        f"The person who told {recipient}, as they stood on {location}, "
        f"the following words was {speaker}."
    )
    
    print(answer_statement)
    print("\nThe quote was:")
    print(quote)

# Execute the function to print the answer
find_quote_origin()