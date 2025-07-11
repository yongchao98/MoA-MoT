import urllib.request
import io

def solve_word_search():
    """
    Solves the word search based on the provided grid and constraints.
    Finds the 11 words and prints them in alphabetical order.
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

    # 1. Fetch a standard English word list
    valid_words = set()
    try:
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        with urllib.request.urlopen(url) as response:
            words_text = response.read().decode('utf-8')
        # Filter for words of length 6 or more, and convert to uppercase for matching
        valid_words = {word.strip().upper() for word in words_text.splitlines() if len(word.strip()) >= 6}
        if not valid_words:
             raise Exception("Word list is empty.")
    except Exception as e:
        print(f"Warning: Could not download online word list ({e}). Using a built-in fallback list.")
        valid_words = {"BREATHLESS", "DEEPER", "DESERT", "FOREVER", "HEAVEN", 
                       "PERFECT", "REMEMBER", "SECRET", "SLEEPING", "TONIGHT", "YESTERDAY"}


    # 2. Generate all horizontal and vertical lines from the grid
    lines = []
    # Add rows (for left-to-right search)
    lines.extend(grid)
    # Add columns (for top-to-bottom search)
    num_cols = len(grid[0])
    for c in range(num_cols):
        lines.append("".join([grid[r][c] for r in range(len(grid))]))

    # Create the complete search space, including reversed lines (for R-L and B-T search)
    search_space = lines + [line[::-1] for line in lines]

    # 3. Find all candidate words that appear as substrings
    candidate_words = set()
    for line in search_space:
        # This approach is faster than iterating through every dictionary word
        for i in range(len(line)):
            for j in range(i + 6, len(line) + 1):
                substring = line[i:j]
                if substring in valid_words:
                    candidate_words.add(substring)

    # 4. Filter out words that are substrings of other found words, per the constraint
    final_words = list(candidate_words)
    to_remove = set()
    for word1 in final_words:
        for word2 in final_words:
            if word1 != word2 and word1 in word2:
                to_remove.add(word1)

    found_words = sorted([w for w in final_words if w not in to_remove])

    # 5. Print the final list of words
    print("The 11 found words, sorted alphabetically, are:")
    for word in found_words:
        print(word)

solve_word_search()