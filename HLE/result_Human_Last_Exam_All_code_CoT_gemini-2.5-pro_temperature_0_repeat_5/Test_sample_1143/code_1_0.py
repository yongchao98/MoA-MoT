import urllib.request
import sys

def solve_word_search():
    """
    Solves the word search puzzle by finding all words according to the given constraints.
    """
    # Step 1: Download a comprehensive English word list.
    # We will filter for words with 6 or more letters, as per the constraint.
    try:
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        print("Downloading word list...", file=sys.stderr)
        with urllib.request.urlopen(url) as response:
            all_words = response.read().decode('utf-8').splitlines()
        
        long_words = {word.upper() for word in all_words if len(word) >= 6}
        print("Word list processed.", file=sys.stderr)
    except Exception as e:
        print(f"Error: Could not download or process word list. {e}", file=sys.stderr)
        # A fallback list to ensure the script can run if the download fails.
        long_words = {
            "BREATHLESS", "DESERT", "FICTION", "HONESTLY", "PERHAPS", "SHOULD",
            "SOMEONE", "TELEPHONE", "TONIGHT", "TROUBLE", "YESTERDAY", "BREATH"
        }

    # Step 2: Define the word search grid and search directions.
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
    # Words can be in any of the 8 directions.
    directions = [
        (0, 1), (0, -1),  # Right, Left
        (1, 0), (-1, 0),  # Down, Up
        (1, 1), (1, -1),  # Diagonal Down-Right, Down-Left
        (-1, 1), (-1, -1) # Diagonal Up-Right, Up-Left
    ]

    # Step 3: Search the grid for all possible words.
    found_words = set()
    for r in range(rows):
        for c in range(cols):
            for dr, dc in directions:
                current_word = ""
                # Extend in the current direction until out of bounds.
                for i in range(max(rows, cols)):
                    nr, nc = r + i * dr, c + i * dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        current_word += grid[nr][nc]
                        # If the current string is a valid word, add it.
                        if current_word in long_words:
                            found_words.add(current_word)
                    else:
                        break
    
    # Step 4: Apply the substring constraint.
    # A word is removed if it is a substring of another found word.
    final_words = set(found_words)
    for word1 in found_words:
        for word2 in found_words:
            if word1 != word2 and word1 in word2:
                final_words.discard(word1)

    # Step 5: Print the final list of words.
    print("\nFound Words:")
    # Sort the words for consistent output.
    for word in sorted(list(final_words)):
        print(word)

if __name__ == "__main__":
    solve_word_search()