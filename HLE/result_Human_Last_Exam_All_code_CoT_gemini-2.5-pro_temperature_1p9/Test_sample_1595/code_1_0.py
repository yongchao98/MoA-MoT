import textwrap

def find_speaker():
    """
    This script identifies and provides context for a quote attributed to a companion of G.K. Chesterton.
    """
    
    # The person who made the statement is identified by scholars.
    speaker = "Maurice Baring"
    
    # Contextual information about the quote.
    recipient = "G.K. Chesterton"
    location = "on the Mount of Olives with a view of Gethsemane in Jerusalem"
    quote = "'Well, anyhow, it must be obvious to anybody that the doctrine of the Fall is the only cheerful view of human life.'"
    source_book = "The New Jerusalem"

    # Explanation text.
    explanation = (
        f"The individual who told {recipient}, while they stood {location}, {quote} was {speaker}.\n\n"
        f"G.K. Chesterton recounts this anecdote in his book '{source_book}'. "
        "While he describes the speaker anonymously as a 'very brilliant and philosophical Anglo-Indian official,' "
        f"it is widely accepted by biographers and scholars that he was referring to his close friend, the diplomat and author, {speaker}."
    )

    # Print the formatted answer.
    print(textwrap.fill(explanation, width=80))

if __name__ == "__main__":
    find_speaker()