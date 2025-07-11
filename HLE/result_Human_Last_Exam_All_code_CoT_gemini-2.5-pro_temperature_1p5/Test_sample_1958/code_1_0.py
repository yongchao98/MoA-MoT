import collections

def solve_bible_in_shakespeare():
    """
    Finds the first book of the Bible, in canonical order, mentioned in a
    Shakespeare play and identifies the play.
    """
    # A list of the Old Testament books in canonical order.
    canonical_books = [
        "Genesis", "Exodus", "Leviticus", "Numbers", "Deuteronomy", "Joshua",
        "Judges", "Ruth", "1 Samuel", "2 Samuel", "1 Kings", "2 Kings",
        "1 Chronicles", "2 Chronicles", "Ezra", "Nehemiah", "Esther", "Job",
        "Psalms", "Proverbs", "Ecclesiastes", "Song of Solomon", "Isaiah",
        "Jeremiah", "Lamentations", "Ezekiel", "Daniel", "Hosea", "Joel",
        "Amos", "Obadiah", "Jonah", "Micah", "Nahum", "Habakkuk",
        "Zephaniah", "Haggai", "Zechariah", "Malachi"
    ]

    # A pre-compiled dictionary of known mentions.
    # The key is the common name of the book, and the value is the play.
    # Note: References to characters (Job, Daniel) are included as the
    # book is named after them. Explicit mentions ("book of Psalms") are also included.
    mentions = collections.OrderedDict([
        ("job", "The Merry Wives of Windsor"),
        # Also in "Henry IV, Part 2", but Merry Wives is a valid instance.
        ("psalms", "Henry IV, Part 2"),
        ("daniel", "The Merchant of Venice"),
        # Other potential mentions are generally ambiguous (e.g., "Kings", "Numbers").
    ])

    # Iterate through the canonical list to find the first match.
    for book in canonical_books:
        # Normalize the book name for searching in our dictionary (e.g., "1 Samuel" -> "samuel")
        search_key = book.lower().split(" ")[-1]

        if search_key in mentions:
            found_book = book
            found_play = mentions[search_key]

            print(f"The first book in the Bible, in canonical order, to be mentioned by name in a Shakespeare play is:")
            print(f"Book: {found_book}")
            print(f"Play: {found_play}")
            return

    print("No biblical book mentions found in the dataset.")

solve_bible_in_shakespeare()