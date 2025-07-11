import requests
import sys

def find_album_from_word_search():
    """
    This script solves the word search puzzle by performing the following steps:
    1. Downloads a standard English dictionary.
    2. Searches the specified grid for words of 6 or more letters in all 4 directions (horizontal, vertical, and their reverses).
    3. Filters the found words to exclude any word that is a substring of another found word.
    4. Prints the final list of 11 words, sorted alphabetically. This list is the key to identifying the album.
    """
    grid = [
        "DESERTGFSG",
        "EEHAJWNLPS",
        "ILONSHIAED",
        "FKUWZEKMEU",
        "ICLHONNECO",
        "RIDKQEISHL",
        "TFIDMVHCLC",
        "EYADRETSEY",
        "PERLBRFPST",
        "BREATHLESS",
    ]
    rows = len(grid)
    cols = len(grid[0])

    try:
        # Using a standard, comprehensive English word list from an online source.
        response = requests.get("https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt")
        response.raise_for_status()
        all_words = response.text.splitlines()
        # Create a set of uppercase words (length >= 6) for fast lookups.
        word_set = {word.upper() for word in all_words if len(word) >= 6}
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the word list. Please check your internet connection. Details: {e}", file=sys.stderr)
        return

    found_words_raw = set()

    # Search horizontally (left-to-right and right-to-left)
    for r in range(rows):
        row_str = grid[r]
        for c in range(cols):
            for length in range(6, cols - c + 1):
                # Check forward substring
                substring = row_str[c:c+length]
                if substring in word_set:
                    found_words_raw.add(substring)
                # Check reverse of the substring
                rev_substring = substring[::-1]
                if rev_substring in word_set:
                    found_words_raw.add(rev_substring)

    # Search vertically (top-to-bottom and bottom-to-top)
    for c in range(cols):
        col_str = "".join(grid[r][c] for r in range(rows))
        for r in range(rows):
            for length in range(6, rows - r + 1):
                # Check forward substring
                substring = col_str[r:r+length]
                if substring in word_set:
                    found_words_raw.add(substring)
                # Check reverse of the substring
                rev_substring = substring[::-1]
                if rev_substring in word_set:
                    found_words_raw.add(rev_substring)

    # Filter out words that are substrings of other found words
    candidates = list(found_words_raw)
    final_words = []
    for word1 in candidates:
        is_substring = False
        for word2 in candidates:
            if word1 != word2 and word1 in word2:
                is_substring = True
                break
        if not is_substring:
            final_words.append(word1)

    # Sort the final list for consistent ordering
    final_words.sort()
    
    if len(final_words) == 11:
        print("Successfully found the 11 words:")
        for word in final_words:
            print(word)
    else:
        print(f"Expected 11 words, but found {len(final_words)}. The grid or dictionary might not match the puzzle's intent.")
        print("Found words:", final_words)

if __name__ == '__main__':
    find_album_from_word_search()