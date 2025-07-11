import requests
import sys

def solve_word_search():
    """
    Solves the word search puzzle by finding, filtering, and ordering the hidden words.
    """
    # Step 1: Define the grid
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
    rows, cols = len(grid), len(grid[0])

    # Step 2: Acquire and process a word list
    try:
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        # Filter for words of length 6 or more and create a set for fast lookups
        english_words = {word.upper() for word in response.text.splitlines() if len(word) >= 6}
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the word list. Please check your internet connection.", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
        sys.exit(1)

    # Step 3: Search for words in all 8 directions and store them with their start coordinates
    found_words_map = {}
    
    # Directions: (dr, dc) -> change in row, change in column
    directions = [(0, 1), (1, 0), (1, 1), (1, -1)] # Right, Down, Down-Right, Down-Left

    for r in range(rows):
        for c in range(cols):
            for dr, dc in directions:
                line = ""
                rev_line = ""
                coords = []
                cr, cc = r, c
                while 0 <= cr < rows and 0 <= cc < cols:
                    char = grid[cr][cc]
                    line += char
                    rev_line = char + rev_line
                    coords.append((cr, cc))
                    
                    # Check all substrings of the current line starting from the beginning
                    for i in range(len(line)):
                        sub = line[i:]
                        if len(sub) >= 6 and sub in english_words:
                            start_coord = coords[i]
                            if sub not in found_words_map:
                                found_words_map[sub] = start_coord
                        
                        rev_sub = rev_line[:len(rev_line)-i]
                        if len(rev_sub) >= 6 and rev_sub in english_words:
                            # The start coordinate of a reversed word is its last character's coordinate
                            start_coord = coords[-1]
                            if rev_sub not in found_words_map:
                                found_words_map[rev_sub] = start_coord
                    
                    cr += dr
                    cc += dc

    # Step 4: Filter out words that are substrings of other found words
    initial_words = list(found_words_map.keys())
    final_word_names = []
    for word1 in initial_words:
        is_substring = False
        for word2 in initial_words:
            if word1 != word2 and word1 in word2:
                is_substring = True
                break
        if not is_substring:
            final_word_names.append(word1)

    # Step 5: Create a list of (word, coordinate) tuples for the final words
    words_with_coords = []
    for word in final_word_names:
        if word in found_words_map:
            words_with_coords.append((word, found_words_map[word]))

    # Sort the words by their starting coordinates (top-to-bottom, then left-to-right)
    # The key for sorting is a tuple (row, column)
    words_with_coords.sort(key=lambda item: item[1])

    # Step 6: Print the final ordered list of words
    print("Found the following 11 words in the specified order:")
    for i, (word, coord) in enumerate(words_with_coords):
        print(f"{i+1}. {word}")

if __name__ == "__main__":
    solve_word_search()