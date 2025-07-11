import requests
import re
import sys

def find_first_bible_book_mention():
    """
    Finds the first book of the Bible, in canonical order, mentioned by name
    in a Shakespeare play.
    """
    # Step 1: Create a list of Bible books in canonical order.
    # We use base names (e.g., "Samuel" for "1 Samuel" and "2 Samuel")
    # as the likely phrasing would be "the book of Samuel".
    bible_books = [
        # Old Testament
        "Genesis", "Exodus", "Leviticus", "Numbers", "Deuteronomy", "Joshua", "Judges", "Ruth",
        "Samuel", "Kings", "Chronicles", "Ezra", "Nehemiah", "Esther", "Job", "Psalms",
        "Proverbs", "Ecclesiastes", "Song of Solomon", "Isaiah", "Jeremiah", "Lamentations",
        "Ezekiel", "Daniel", "Hosea", "Joel", "Amos", "Obadiah", "Jonah", "Micah", "Nahum",
        "Habakkuk", "Zephaniah", "Haggai", "Zechariah", "Malachi",
        # New Testament
        "Matthew", "Mark", "Luke", "John", "Acts", "Romans", "Corinthians", "Galatians",
        "Ephesians", "Philippians", "Colossians", "Thessalonians", "Timothy", "Titus",
        "Philemon", "Hebrews", "James", "Peter", "Jude", "Revelation"
    ]

    # Step 2: Download Shakespeare's complete works.
    try:
        url = "https://www.gutenberg.org/cache/epub/100/pg100.txt"
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        # Use 'utf-8-sig' to handle the potential BOM (Byte Order Mark) at the start of the file
        shakespeare_text = response.content.decode('utf-8-sig')
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the text file. {e}", file=sys.stderr)
        return

    # Step 3: Split the text into individual plays.
    # This list of titles helps to reliably separate the plays from the Gutenberg text.
    play_titles = [
        "THE TEMPEST", "THE TWO GENTLEMEN OF VERONA", "THE MERRY WIVES OF WINDSOR",
        "MEASURE FOR MEASURE", "THE COMEDY OF ERRORS", "MUCH ADO ABOUT NOTHING",
        "LOVE'S LABOUR'S LOST", "A MIDSUMMER NIGHT'S DREAM", "THE MERCHANT OF VENICE",
        "AS YOU LIKE IT", "THE TAMING OF THE SHREW", "ALL'S WELL THAT ENDS WELL",
        "TWELFTH NIGHT; OR, WHAT YOU WILL", "THE WINTER'S TALE", "THE LIFE AND DEATH OF KING JOHN",
        "THE TRAGEDY OF KING RICHARD THE SECOND", "THE FIRST PART OF KING HENRY THE FOURTH",
        "THE SECOND PART OF KING HENRY THE FOURTH", "THE LIFE OF KING HENRY THE FIFTH",
        "THE FIRST PART OF KING HENRY THE SIXTH", "THE SECOND PART OF KING HENRY THE SIXTH",
        "THE THIRD PART OF KING HENRY THE SIXTH", "THE TRAGEDY OF KING RICHARD THE THIRD",
        "THE FAMOUS HISTORY OF THE LIFE OF KING HENRY THE EIGHTH", "TROILUS AND CRESSIDA",
        "THE TRAGEDY OF CORIOLANUS", "TITUS ANDRONICUS", "THE TRAGEDY OF ROMEO AND JULIET",
        "TIMON OF ATHENS", "THE TRAGEDY OF JULIUS CAESAR", "THE TRAGEDY OF MACBETH",
        "THE TRAGEDY OF HAMLET, PRINCE OF DENMARK", "KING LEAR", "OTHELLO, THE MOOR OF VENICE",
        "ANTONY AND CLEOPATRA", "CYMBELINE"
    ]

    plays = {}
    current_pos = 0
    for i, title in enumerate(play_titles):
        start_index = shakespeare_text.find(title, current_pos)
        if start_index == -1:
            continue

        end_index = -1
        if i + 1 < len(play_titles):
            # Find the start of the next play to determine the end of the current one
            next_title = play_titles[i+1]
            end_index = shakespeare_text.find(next_title, start_index + len(title))

        if end_index != -1:
            plays[title] = shakespeare_text[start_index:end_index]
        else: # Handle the last play in the list
            plays[title] = shakespeare_text[start_index:]
        current_pos = start_index + 1


    # Step 4: Iterate through books and plays to find the first mention.
    for book in bible_books:
        for play_title, play_text in plays.items():
            # Search for the pattern "the book of [Book Name]" case-insensitively.
            # \b ensures we match whole words only.
            search_pattern = re.compile(fr'\bthe book of {book}\b', re.IGNORECASE)
            match = search_pattern.search(play_text)

            if match:
                # If a match is found, we have our answer.
                print(f"The first book of the Bible mentioned by name, in canonical order, is the Book of {book}.")
                # .title() makes the play title look nice, e.g., "The Merry Wives Of Windsor"
                print(f"It is mentioned in Shakespeare's play: '{play_title.title()}'.")
                return # Exit after finding the first one.

    print("No explicit mention of a Bible book by name in the format 'the book of...' was found.")

if __name__ == '__main__':
    find_first_bible_book_mention()