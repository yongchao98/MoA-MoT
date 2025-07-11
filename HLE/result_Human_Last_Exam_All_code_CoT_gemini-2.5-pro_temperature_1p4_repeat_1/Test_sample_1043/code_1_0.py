def find_poet_and_emphasized_word():
    """
    Identifies the poet of the provided lines and the word they emphasized
    most in a spoken performance.

    The lines are from Frank O'Hara's poem "This is a Portrait of V.R. Lang".
    In recordings of O'Hara reading this poem, the word "holler" is delivered
    with a distinct and notable emphasis, making it stand out from the
    surrounding words.
    """
    poet = "O'Hara"
    emphasized_word = "holler"
    
    # Format the answer as "LastName, word"
    formatted_answer = f"{poet}, {emphasized_word}"
    
    print(formatted_answer)

find_poet_and_emphasized_word()