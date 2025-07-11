def find_speaker():
    """
    This function provides the answer to the user's question based on G.K. Chesterton's writings.
    The person who made the statement to Chesterton is not named in the source text.
    Instead, Chesterton provides a description.
    """
    
    # Information gathered from G.K. Chesterton's book "The Everlasting Man", Part I, Chapter 5.
    description = "an old and distinguished diplomatist and man of the world"
    source = "G.K. Chesterton's book, 'The Everlasting Man'"

    # Format and print the answer
    output = (
        f"The person who said this to G.K. Chesterton is not explicitly named in his writings.\n"
        f"In the source of the anecdote, {source}, Chesterton describes the speaker as:\n"
        f"'{description}'."
    )
    
    print(output)

find_speaker()