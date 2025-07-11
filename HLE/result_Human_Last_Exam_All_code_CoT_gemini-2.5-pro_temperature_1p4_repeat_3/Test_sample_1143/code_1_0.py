import os

def solve_word_search():
    """
    Finds 11 specific words in a grid based on given rules.
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
        "BREATHLESS"
    ]

    # --- Step 1: Load a dictionary ---
    # Using a common dictionary path. Users on other systems might need to change this path.
    # A comprehensive word list is needed for best results.
    dictionary_path = '/usr/share/dict/words'
    if not os.path.exists(dictionary_path):
        print(f"Dictionary file not found at {dictionary_path}.")
        print("Please modify the 'dictionary_path' variable to point to a word list file.")
        return
        
    with open(dictionary_path, 'r') as f:
        dictionary = {line.strip().lower() for line in f}

    # --- Step 2: Generate all possible lines from the grid ---
    lines = []
    # Rows (L->R)
    lines.extend(grid)
    # Rows (R->L)
    lines.extend([row[::-1] for row in grid])
    # Columns (T->B)
    cols = [''.join([grid[i][j] for i in range(len(grid))]) for j in range(len(grid[0]))]
    lines.extend(cols)
    # Columns (B->T)
    lines.extend([col[::-1] for col in cols])
    
    # --- Step 3: Find all possible valid words ---
    potential_words = set()
    # To maintain a deterministic "found" order, we use a list first.
    ordered_potential_words = []

    for line in lines:
        for i in range(len(line)):
            for j in range(i + 6, len(line) + 1):
                word = line[i:j].lower()
                if word in dictionary:
                    if word not in potential_words:
                        ordered_potential_words.append(word)
                        potential_words.add(word)

    # --- Step 4: Filter out words that are substrings of other found words ---
    words_to_remove = set()
    for word1 in potential_words:
        for word2 in potential_words:
            if word1 != word2 and word1 in word2:
                words_to_remove.add(word1)
    
    final_words = [word for word in ordered_potential_words if word not in words_to_remove]
    
    # --- Step 5: Print the final list of words ---
    print("The 11 found words in order of discovery are:")
    for word in final_words:
        print(word.upper())

solve_word_search()
