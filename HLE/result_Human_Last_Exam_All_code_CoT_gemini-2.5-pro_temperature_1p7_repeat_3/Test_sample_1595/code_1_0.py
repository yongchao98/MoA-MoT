def find_speaker():
    """
    This function identifies and prints the name of the person
    who spoke a specific quote to G.K. Chesterton.
    """
    
    # The known information based on historical and literary sources.
    speaker = "his wife, Frances Blogg Chesterton"
    recipient = "G.K. Chesterton"
    location = "the Mount of Olives"
    quote = "Well, anyhow, it must be obvious to anybody that the doctrine of the Fall is the only cheerful view of human life."

    # Print the answer to the user's question.
    print(f"It was {speaker} who told {recipient}, as they stood on {location}, the following:")
    print(f"'{quote}'")

# Execute the function to provide the answer.
find_speaker()