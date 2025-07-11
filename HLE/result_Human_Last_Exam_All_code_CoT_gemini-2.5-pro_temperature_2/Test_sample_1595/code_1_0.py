def find_the_speaker():
    """
    This function identifies and prints the person who spoke the quote
    to G.K. Chesterton on the Mount of Olives. The information is found
    in Chesterton's own writing.
    """
    speaker = "his wife, Frances Blogg Chesterton"
    listener = "G.K. Chesterton"
    location = "the Mount of Olives with a view of Gethsemane in Jerusalem"
    quote = "'Well, anyhow, it must be obvious to anybody that the doctrine of the Fall is the only cheerful view of human life'"

    # Construct the final sentence providing the answer.
    answer_text = (
        f"According to G.K. Chesterton's own account in his book 'The Everlasting Man', "
        f"it was {speaker} who told him, as they stood on {location}, "
        f"{quote}."
    )
    
    print(answer_text)

find_the_speaker()