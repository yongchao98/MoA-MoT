def find_first_bible_book_mention():
    """
    Finds the first book of the Bible, in canonical order,
    that is mentioned by name in a Shakespeare play.
    """

    # A list of Old Testament books in canonical order. We only need to go
    # as far as the first confirmed mention.
    bible_books_in_order = [
        "Genesis", "Exodus", "Leviticus", "Numbers", "Deuteronomy",
        "Joshua", "Judges", "Ruth", "1 Samuel", "2 Samuel", "1 Kings",
        "2 Kings", "1 Chronicles", "2 Chronicles", "Ezra", "Nehemiah",
        "Esther", "Job", "Psalms", "Proverbs"
    ]

    # A dictionary of confirmed, explicit mentions of Biblical books by name
    # in Shakespeare's plays. The name "Job" is used, which is the title of the
    # book. This is the earliest book in the canon to be named.
    # The reference to "Psalms" is also valid but occurs later in the canon.
    mentions = {
        "Job": {
            "play": "Henry IV, Part 2",
            "quote": "I am as poor as Job, my lord, but not so patient.",
            "location": "Act I, Scene 2"
        },
        "Psalms": {
            "play": "The Merry Wives of Windsor",
            "quote": "they do no more adhere and keep pace together than the Hundredth Psalm to the tune of ‘Green Sleeves.’",
            "location": "Act II, Scene 1"
        }
    }

    # Iterate through the canonical list and stop at the first book found in our mentions
    for book in bible_books_in_order:
        if book in mentions:
            found_book = book
            found_play = mentions[book]["play"]
            
            print(f"The first book of the Bible, in canonical order, mentioned by name in a Shakespeare play is the Book of {found_book}.")
            print(f"It is mentioned in the play '{found_play}'.")
            return

find_first_bible_book_mention()