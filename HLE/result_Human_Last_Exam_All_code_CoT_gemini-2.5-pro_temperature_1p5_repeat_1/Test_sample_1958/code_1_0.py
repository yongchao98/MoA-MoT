import nltk
import re
import ssl

def find_first_bible_book_in_shakespeare():
    """
    Finds the first book of the Bible, in canonical order, mentioned by name
    in a Shakespeare play.
    """
    # --- Step 1: Prepare NLTK Data ---
    # Handle potential SSL errors when downloading NLTK data
    try:
        _create_unverified_https_context = ssl._create_unverified_context
    except AttributeError:
        pass
    else:
        ssl._create_default_https_context = _create_unverified_https_context

    try:
        nltk.data.find('corpora/gutenberg')
    except LookupError:
        print("Downloading 'gutenberg' corpus from NLTK...")
        nltk.download('gutenberg', quiet=True)
        print("Download complete.")

    # --- Step 2: Define Bible Books and Search Logic ---
    # List of Bible books in canonical order
    bible_books = [
        "Genesis", "Exodus", "Leviticus", "Numbers", "Deuteronomy", "Joshua",
        "Judges", "Ruth", "Samuel", "Kings", "Chronicles", "Ezra", "Nehemiah",
        "Esther", "Job", "Psalms", "Proverbs", "Ecclesiastes", "Song of Solomon",
        "Isaiah", "Jeremiah", "Lamentations", "Ezekiel", "Daniel", "Hosea",
        "Joel", "Amos", "Obadiah", "Jonah", "Micah", "Nahum", "Habakkuk",
        "Zephaniah", "Haggai", "Zechariah", "Malachi", "Matthew", "Mark", "Luke",
        "John", "Acts", "Romans", "Corinthians", "Galatians", "Ephesians",
        "Philippians", "Colossians", "Thessalonians", "Timothy", "Titus",
        "Philemon", "Hebrews", "James", "Peter", "Jude", "Revelation"
    ]

    # Books whose names are common words, requiring a more specific search
    common_word_books = {
        "Numbers", "Judges", "Kings", "Chronicles", "Proverbs", "Lamentations",
        "Mark", "Luke", "John", "Acts", "James", "Jude"
    }

    # --- Step 3: Search Through Shakespeare's Plays ---
    shakespeare_files = [f for f in nltk.corpus.gutenberg.fileids() if f.startswith('shakespeare-')]

    found_book = None
    found_play = None
    
    for book in bible_books:
        # Create a search pattern based on whether the book name is a common word
        if book in common_word_books:
            # Search for the more explicit "book of..." phrase
            pattern = re.compile(r'book of ' + re.escape(book), re.IGNORECASE)
        else:
            # For less common names, search for the name as a whole word
            pattern = re.compile(r'\b' + re.escape(book) + r'\b', re.IGNORECASE)

        for play_file in shakespeare_files:
            text = nltk.corpus.gutenberg.raw(play_file)
            if pattern.search(text):
                found_book = book
                # Format play name for readability
                play_name = play_file.replace('shakespeare-', '').replace('.txt', '').replace('-', ' ').title()
                found_play = play_name
                break # Stop searching plays for this book
        
        if found_book:
            break # Stop searching books

    # --- Step 4: Print the Result ---
    if found_book and found_play:
        print(f"The first book in the Bible, in canonical order, to be mentioned by name in a Shakespeare play is '{found_book}'.")
        print(f"It is mentioned in the play: '{found_play}'.")
    else:
        print("Could not find any biblical book references in Shakespeare's plays with this method.")

if __name__ == '__main__':
    find_first_bible_book_in_shakespeare()