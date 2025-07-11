import requests
import sys

def find_album_from_word_search():
    """
    This script solves the word search puzzle to find 11 specific English words.
    These words are the key to identifying the secret album.
    """
    
    # Step 1: Define the grid and search directions
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
    directions = [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]

    # Step 2: Fetch a comprehensive English word list for validation
    word_set = set()
    try:
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        word_set = set(word.strip().lower() for word in response.text.splitlines())
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not fetch the required English word list from the internet: {e}", file=sys.stderr)
        print("Please check your internet connection and try again.", file=sys.stderr)
        return

    # Step 3: Find all possible words in the grid that are 6+ letters long
    all_possible_words = set()
    for r in range(rows):
        for c in range(cols):
            for dr, dc in directions:
                word = ""
                current_r, current_c = r, c
                while 0 <= current_r < rows and 0 <= current_c < cols:
                    word += grid[current_r][current_c]
                    if len(word) >= 6:
                        if word.lower() in word_set:
                            all_possible_words.add(word.lower())
                    current_r += dr
                    current_c += dc

    # Step 4: Filter the words based on the constraint that no word is a substring
    # of a longer found word.
    final_words = set()
    for word1 in all_possible_words:
        is_substring_of_longer = False
        for word2 in all_possible_words:
            if len(word1) < len(word2) and word1 in word2:
                is_substring_of_longer = True
                break
        if not is_substring_of_longer:
            final_words.add(word1)
            
    # Step 5: Print the final list of 11 found words, which are the clues to the album.
    # The prompt mentioned an "equation," which is likely a template error.
    # Printing the words themselves is the correct interpretation.
    sorted_words = sorted(list(final_words))
    
    print("The 11 found words are:")
    for word in sorted_words:
        print(word)

if __name__ == '__main__':
    find_album_from_word_search()
