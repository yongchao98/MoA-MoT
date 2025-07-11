def find_speaker():
    """
    This function identifies and prints the name of the person who spoke a specific
    quote to G.K. Chesterton.
    """
    # The quote in question
    quote = "Well, anyhow, it must be obvious to anybody that the doctrine of the Fall is the only cheerful view of human life."

    # The person who received the quote
    recipient = "G.K. Chesterton"

    # The location where the quote was said
    location = "the Mount of Olives with a view of Gethsemane in Jerusalem"

    # According to G.K. Chesterton's "Autobiography" (Chapter 14), the speaker was his wife.
    speaker = "Frances Chesterton (his wife)"

    # Print the answer
    print(f"The person who told {recipient}, while on {location}, '{quote}' was:")
    print(speaker)

if __name__ == "__main__":
    find_speaker()