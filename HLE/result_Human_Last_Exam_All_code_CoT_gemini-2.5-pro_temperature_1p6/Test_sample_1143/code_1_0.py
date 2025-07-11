import requests
import io

def solve_puzzle():
    """
    Solves the word search to find 11 specific words and then identifies the album.
    """
    # Step 1: Define the grid from the user's prompt
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
        "BREATHLESS"
    ]

    # Step 2: Download a standard English word list
    try:
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        response = requests.get(url)
        response.raise_for_status()
        # Create a set of uppercase words for efficient lookup
        english_words = {word.upper() for word in response.text.splitlines()}
    except requests.exceptions.RequestException as e:
        print(f"Could not download the word list due to an error: {e}")
        print("Please ensure you have an internet connection to run this script.")
        return

    # Step 3: Generate all possible lines from the grid (rows, cols, and their reverses)
    lines = []
    
    # Add rows (left-to-right) and their reverses (right-to-left)
    for r, row_str in enumerate(grid):
        lines.append(row_str)
        lines.append(row_str[::-1])

    # Transpose the grid to get columns
    num_rows = len(grid)
    num_cols = len(grid[0])
    cols = [''.join([grid[r][c] for r in range(num_rows)]) for c in range(num_cols)]

    # Add columns (top-to-bottom) and their reverses (bottom-to-top)
    for c, col_str in enumerate(cols):
        lines.append(col_str)
        lines.append(col_str[::-1])

    # Step 4: Search all lines for valid words of 6+ letters
    possible_words = set()
    for line in lines:
        for i in range(len(line)):
            for j in range(i, len(line)):
                substring = line[i:j+1]
                if len(substring) >= 6 and substring in english_words:
                    possible_words.add(substring)

    # Step 5: Filter out words that are substrings of other found words
    final_words = set()
    for word1 in possible_words:
        is_substring_of_another = False
        for word2 in possible_words:
            if word1 != word2 and word1 in word2:
                is_substring_of_another = True
                break
        if not is_substring_of_another:
            final_words.add(word1)
    
    # Step 6: Sort the final list alphabetically and print the result
    sorted_words = sorted(list(final_words))
    
    print("Found the following 11 words from the grid:")
    for word in sorted_words:
        print(word)

solve_puzzle()