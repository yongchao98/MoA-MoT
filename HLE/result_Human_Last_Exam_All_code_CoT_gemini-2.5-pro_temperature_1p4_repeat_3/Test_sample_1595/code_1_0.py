def find_speaker():
    """
    This function identifies and prints the name of the person who spoke a specific
    quote to G.K. Chesterton.
    """
    # The person who spoke the words was G.K. Chesterton's brother.
    speaker = "Cecil Chesterton"

    # The person being spoken to.
    listener = "G.K. Chesterton"

    # The location of the event.
    location = "the Mount of Olives with a view of Gethsemane in Jerusalem"

    # The quote.
    quote = "'Well, anyhow, it must be obvious to anybody that the doctrine of the Fall is the only cheerful view of human life.'"

    # Print the final answer.
    print(f"It was {speaker} who told {listener}, as they both stood on {location}, {quote}.")


if __name__ == "__main__":
    find_speaker()