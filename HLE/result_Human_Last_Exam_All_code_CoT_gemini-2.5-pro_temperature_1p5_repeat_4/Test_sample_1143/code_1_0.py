import requests

def solve_word_search_puzzle():
    """
    Finds words in the provided grid according to the puzzle's rules.
    """
    # The word search grid provided by the user.
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

    print("Starting word search on the provided grid...")

    # Step 1: Fetch a comprehensive list of English words.
    try:
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        # Create a set of uppercase words (length >= 6) for efficient lookup.
        english_words = {word.strip().upper() for word in response.text.splitlines() if len(word.strip()) >= 6}
        print(f"Successfully loaded {len(english_words)} English words.")
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not fetch the word list from the internet. {e}")
        return

    # Step 2: Generate all possible lines to search from the grid (rows, cols, and their reverses).
    lines_to_search = []
    # Horizontal lines
    for row_str in grid:
        lines_to_search.append(row_str)
        lines_to_search.append(row_str[::-1])
    # Vertical lines
    num_rows = len(grid)
    num_cols = len(grid[0])
    for c in range(num_cols):
        col_str = "".join([grid[r][c] for r in range(num_rows)])
        lines_to_search.append(col_str)
        lines_to_search.append(col_str[::-1])

    # Step 3: Find all substrings in those lines that are valid English words.
    candidate_words = set()
    for line in lines_to_search:
        for i in range(len(line)):
            for j in range(i + 6, len(line) + 1):
                substring = line[i:j]
                if substring in english_words:
                    candidate_words.add(substring)

    # Step 4: Filter out words that are substrings of other longer found words.
    final_words = set(candidate_words)
    for word1 in candidate_words:
        for word2 in candidate_words:
            if word1 != word2 and word1 in word2:
                final_words.discard(word1)
                break
    
    # Step 5: Print the final list of words found.
    print("\n--- Found Words ---")
    if not final_words:
        print("No words matching the criteria were found.")
    else:
        # Sort for consistent output
        sorted_words = sorted(list(final_words))
        for word in sorted_words:
            print(word)
    
    print(f"\nFound a total of {len(final_words)} words.")


solve_word_search_puzzle()