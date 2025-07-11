def find_speaker():
    """
    This function identifies and prints the name of the person
    who spoke a specific quote to G.K. Chesterton in Jerusalem.
    """
    speaker = "H. A. L. Fisher"
    listener = "G.K. Chesterton"
    location = "the Mount of Olives in Jerusalem"
    quote = "Well, anyhow, it must be obvious to anybody that the doctrine of the Fall is the only cheerful view of human life"

    # In his book "The New Jerusalem", Chesterton describes this person as a
    # "distinguished... diplomatist", "Member of Parliament", and "Cabinet Minister".
    # Historians have identified this person as H. A. L. Fisher.

    print(f"The person who told G.K. Chesterton, '{quote}', while on {location} was {speaker}.")

if __name__ == "__main__":
    find_speaker()