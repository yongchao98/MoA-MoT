def solve_chesterton_quote_riddle():
    """
    This script identifies and prints the name of the person who spoke a specific
    quote to G.K. Chesterton in Jerusalem.
    """
    
    # The quote is from G.K. Chesterton's book "The Everlasting Man".
    # The person who said it was his wife.
    speaker_name = "his wife, Frances Blogg Chesterton"
    
    # Constructing the full answer
    answer_text = (
        f"The person who told G.K. Chesterton, 'Well, anyhow, it must be obvious to anybody "
        f"that the doctrine of the Fall is the only cheerful view of human life,' "
        f"was {speaker_name}."
    )
    
    print(answer_text)

if __name__ == "__main__":
    solve_chesterton_quote_riddle()