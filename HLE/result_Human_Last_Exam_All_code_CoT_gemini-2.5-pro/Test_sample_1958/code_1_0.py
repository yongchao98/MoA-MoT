import requests
import re
import ssl

# Bypass SSL verification for environments that might need it.
# This is sometimes necessary for accessing older resources like Gutenberg.
try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    pass
else:
    ssl._create_default_https_context = _create_unverified_https_context

def find_first_bible_book_in_shakespeare():
    """
    Finds the first book of the Bible, in canonical order, mentioned by name
    in a Shakespeare play.
    """

    # List of Shakespeare's plays from Project Gutenberg by ID
    plays = {
        "Henry VI, Part 2": 1113,
        "Henry VI, Part 3": 1114,
        "Henry VI, Part 1": 1112,
        "Richard III": 1503,
        "The Comedy of Errors": 23046,
        "Titus Andronicus": 1789,
        "The Taming of the Shrew": 1508,
        "The Two Gentlemen of Verona": 1111,
        "Love's Labour's Lost": 23044,
        "Romeo and Juliet": 1513,
        "Richard II": 1780,
        "A Midsummer Night's Dream": 1514,
        "King John": 1518,
        "The Merchant of Venice": 2243,
        "Henry IV, Part 1": 1782,
        "Henry IV, Part 2": 1783,
        "Much Ado About Nothing": 2240,
        "Henry V": 2252,
        "Julius Caesar": 2263,
        "As You Like It": 1786,
        "Twelfth Night": 1525,
        "Hamlet": 1524,
        "The Merry Wives of Windsor": 1779,
        "Troilus and Cressida": 1790,
        "Othello": 2267,
        "All's Well That Ends Well": 1526,
        "Measure for Measure": 1121,
        "King Lear": 2265,
        "Macbeth": 2264,
        "Antony and Cleopatra": 2266,
        "Coriolanus": 1794,
        "Timon of Athens": 1793,
        "Pericles, Prince of Tyre": 2246,
        "Cymbeline": 2248,
        "The Winter's Tale": 1539,
        "The Tempest": 23042,
        "Henry VIII": 1117,
    }

    # Canonical order of Bible books
    bible_books = [
        "Genesis", "Exodus", "Leviticus", "Numbers", "Deuteronomy", "Joshua",
        "Judges", "Ruth", "Samuel", "Kings", "Chronicles", "Ezra", "Nehemiah",
        "Esther", "Job", "Psalms", "Proverbs", "Ecclesiastes", "Song of Solomon",
        "Isaiah", "Jeremiah", "Lamentations", "Ezekiel", "Daniel", "Hosea",
        "Joel", "Amos", "Obadiah", "Jonah", "Micah", "Nahum", "Habakkuk",
        "Zephaniah", "Haggai", "Zechariah", "Malachi", "Matthew", "Mark",
        "Luke", "John", "Acts", "Romans", "Corinthians", "Galatians",
        "Ephesians", "Philippians", "Colossians", "Thessalonians", "Timothy",
        "Titus", "Philemon", "Hebrews", "James", "Peter", "Jude", "Revelation"
    ]

    for book in bible_books:
        # The pattern looks for the book name as a whole word, case-insensitive
        pattern = re.compile(r'\b' + re.escape(book) + r'\b', re.IGNORECASE)

        for play_title, play_id in plays.items():
            try:
                # Fetch play text from Project Gutenberg
                url = f"https://www.gutenberg.org/cache/epub/{play_id}/pg{play_id}.txt"
                response = requests.get(url, timeout=10)
                response.raise_for_status()
                text = response.text
            except requests.exceptions.RequestException as e:
                # Silently continue if a single play fails to download
                continue

            # Split text into lines to find the exact line
            lines = text.splitlines()
            for i, line in enumerate(lines):
                match = pattern.search(line)
                if match:
                    # To provide context, find the speaker if possible
                    speaker = ""
                    # Look backwards for a line in all caps ending with a period, a common speaker format
                    for j in range(i - 1, max(i - 10, -1), -1):
                        prev_line = lines[j].strip()
                        if prev_line.isupper() and prev_line.endswith(('.', ' ')):
                            speaker = prev_line.strip('. ')
                            break
                    
                    print(f"The first book of the Bible mentioned in Shakespeare's plays is '{book}'.")
                    print(f"It is mentioned in the play '{play_title}'.")
                    print("\nHere is the relevant passage:")
                    print("-" * 30)
                    if speaker:
                        print(f"Speaker: {speaker}")
                    print(f"Line: {line.strip()}")
                    print("-" * 30)
                    return

find_first_bible_book_in_shakespeare()