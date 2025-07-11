def find_answer():
    """
    This function provides the answer to the user's trivia question
    about G.K. Chesterton.
    """
    # The person who accompanied Chesterton and made the remark is known from
    # biographical sources.
    speaker = "Maurice Baring"
    
    # The person being spoken to.
    listener = "G.K. Chesterton"
    
    # The quote in question.
    quote = "'Well, anyhow, it must be obvious to anybody that the doctrine of the Fall is the only cheerful view of human life'"

    # We construct and print the final answer.
    print(f"It was {speaker} who told {listener}, as they stood on the Mount of Olives, {quote}.")

find_answer()