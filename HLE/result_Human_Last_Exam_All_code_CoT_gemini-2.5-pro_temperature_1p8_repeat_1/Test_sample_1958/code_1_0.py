import re
import requests

def find_first_bible_book_in_shakespeare():
    """
    Finds the first book of the Bible, in canonical order, mentioned by name
    in Shakespeare's complete works and identifies the play.
    """
    # A list of the books of the Old Testament in canonical order.
    # We only need to go far enough to find the first match.
    BIBLE_BOOKS_OT = [
        "Genesis", "Exodus", "Leviticus", "Numbers", "Deuteronomy",
        "Joshua", "Judges", "Ruth", "Samuel", "Kings", "Chronicles",
        "Ezra", "Nehemiah", "Esther", "Job", "Psalms", "Proverbs",
        "Ecclesiastes", "Song of Solomon", "Isaiah", "Jeremiah",
        "Lamentations", "Ezekiel", "Daniel", "Hosea", "Joel", "Amos",
        "Obadiah", "Jonah", "Micah", "Nahum", "Habakkuk", "Zephaniah",
        "Haggai", "Zechariah", "Malachi"
    ]

    # URL for the plain text of Shakespeare's complete works
    url = "https://www.gutenberg.org/files/100/100-0.txt"
    
    try:
        response = requests.get(url)
        # Raise an HTTPError for bad responses (4xx or 5xx)
        response.raise_for_status()
        text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Failed to download Shakespeare's works. Error: {e}")
        return

    # First, find all potential play titles and their starting positions.
    # Titles in the Gutenberg file are typically all-caps, surrounded by newlines.
    # We filter for titles with more than one word and avoid common non-title headers.
    title_pattern = re.compile(r'\n{2,}\s*([A-Z\s,]{5,})\s*\n{2,}')
    titles_with_indices = []
    for match in title_pattern.finditer(text):
        title = match.group(1).strip()
        # Filter out common headers that are not play titles
        if "ACT" not in title and "SCENE" not in title and "PROLOGUE" not in title \
           and "EPILOGUE" not in title and "CONTENTS" not in title and "FINIS" not in title \
           and len(title.split()) > 1:
            titles_with_indices.append((title, match.start()))

    # Iterate through the books of the Bible in order
    for book in BIBLE_BOOKS_OT:
        # Create a search pattern for the book name (as a whole word, case-insensitive).
        # We add a special case for "Psalms" to also match "Psalm".
        if book == "Psalms":
            search_pattern_str = r'\b(Psalm|Psalms)\b'
        else:
            search_pattern_str = r'\b' + re.escape(book) + r'\b'
        
        book_pattern = re.compile(search_pattern_str, re.IGNORECASE)
        search_result = book_pattern.search(text)
        
        # If a match is found, we have our answer
        if search_result:
            found_book = book
            match_index = search_result.start()
            
            # Now, determine which play this mention belongs to by finding
            # the last title that appeared before the match index.
            current_play_title = "Unknown Play"
            for title, title_index in titles_with_indices:
                if title_index < match_index:
                    current_play_title = title
                else:
                    # We have gone past the match, so the previous title is correct.
                    break
            
            # Format the title nicely for printing (Title Case)
            formatted_title = ' '.join(word.capitalize() for word in current_play_title.lower().split())
            
            print(f"The first book of the Bible, in canonical order, to be mentioned by name in a Shakespeare play is '{found_book}'.")
            print(f"It is mentioned in the play '{formatted_title}'.")
            
            return # Exit after finding the first one

    print("No mentions of Bible books were found in Shakespeare's works.")

if __name__ == '__main__':
    find_first_bible_book_in_shakespeare()