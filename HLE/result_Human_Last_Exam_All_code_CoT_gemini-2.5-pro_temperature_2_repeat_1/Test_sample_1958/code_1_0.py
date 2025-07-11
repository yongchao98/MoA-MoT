import collections

def find_first_mentioned_bible_book():
    """
    Finds the first book of the Bible, in canonical order,
    mentioned by name in a Shakespeare play.
    """
    # Canonical order of the 66 books of the Protestant Bible
    canonical_bible_books = [
        "Genesis", "Exodus", "Leviticus", "Numbers", "Deuteronomy",
        "Joshua", "Judges", "Ruth", "1 Samuel", "2 Samuel", "1 Kings", "2 Kings",
        "1 Chronicles", "2 Chronicles", "Ezra", "Nehemiah", "Esther", "Job",
        "Psalms", "Proverbs", "Ecclesiastes", "Song of Solomon", "Isaiah",
        "Jeremiah", "Lamentations", "Ezekiel", "Daniel", "Hosea", "Joel",
        "Amos", "Obadiah", "Jonah", "Micah", "Nahum", "Habakkuk", "Zephaniah",
        "Haggai", "Zechariah", "Malachi", "Matthew", "Mark", "Luke", "John",
        "Acts", "Romans", "1 Corinthians", "2 Corinthians", "Galatians",
        "Ephesians", "Philippians", "Colossians", "1 Thessalonians",
        "2 Thessalonians", "1 Timothy", "2 Timothy", "Titus", "Philemon",
        "Hebrews", "James", "1 Peter", "2 Peter", "1 John", "2 John",
        "3 John", "Jude", "Revelation"
    ]

    # A dictionary of direct Bible book mentions in Shakespeare's plays.
    # We use a key like "Corinthians" to represent both 1 & 2 Corinthians.
    # The most unambiguous mentions are included.
    mentions_in_shakespeare = {
        "Psalms": "The Merry Wives of Windsor",
        "Corinthians": "Henry IV, Part 1",
        # Note: "Daniel" and "Job" are mentioned as people, not books.
        # "Proverbs" is mentioned, but ambiguously ("patch grief with proverbs").
        # The mention of "the Hundredth Psalm" is a very direct reference to the Book of Psalms.
    }
    
    # We use a special key for Corinthians since it's mentioned without 1 or 2
    mentions_in_shakespeare["1 Corinthians"] = mentions_in_shakespeare["Corinthians"]
    mentions_in_shakespeare["2 Corinthians"] = mentions_in_shakespeare["Corinthians"]

    # Iterate through the canonical list to find the first match
    for book in canonical_bible_books:
        if book in mentions_in_shakespeare:
            play = mentions_in_shakespeare[book]
            # Found the first one
            print(f"The first book in the Bible, in canonical order, mentioned by name in a Shakespeare play is the Book of {book}.")
            print(f"It is mentioned in the play: '{play}'.")
            return

find_first_mentioned_bible_book()