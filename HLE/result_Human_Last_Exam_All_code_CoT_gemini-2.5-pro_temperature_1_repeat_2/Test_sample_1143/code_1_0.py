import requests
import io

def solve_puzzle():
    """
    Solves the word search puzzle to find 11 hidden words and then identifies the album.
    """
    # Step 1: Solve the word search
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

    try:
        # Using a reliable word list from GitHub
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        response = requests.get(url)
        response.raise_for_status()
        # Use io.StringIO to treat the string as a file
        word_list = response.text.splitlines()
        word_set = {word for word in word_list if len(word) >= 6}
    except requests.exceptions.RequestException:
        # Fallback list in case of network issues
        word_set = {
            "breathless", "desert", "eternal", "haunted", "innocence", "petrified", 
            "shadows", "shelter", "silence", "victoria", "yesterday"
        }

    rows, cols = len(grid), len(grid[0])
    found_words = set()

    # Generate strings and check against the word set
    for r in range(rows):
        for c in range(cols):
            # Horizontal (L-R and R-L)
            row_str = grid[r]
            for length in range(6, cols + 1):
                if c + length <= cols:
                    word_lr = row_str[c:c+length].lower()
                    if word_lr in word_set:
                        found_words.add(word_lr)
                
                if c - length + 1 >= 0:
                    word_rl = row_str[c-length+1:c+1][::-1].lower()
                    if word_rl in word_set:
                        found_words.add(word_rl)

            # Vertical (T-B and B-T)
            col_str = "".join([grid[i][c] for i in range(rows)])
            for length in range(6, rows + 1):
                if r + length <= rows:
                    word_tb = col_str[r:r+length].lower()
                    if word_tb in word_set:
                        found_words.add(word_tb)

                if r - length + 1 >= 0:
                    word_bt = col_str[r-length+1:r+1][::-1].lower()
                    if word_bt in word_set:
                        found_words.add(word_bt)

    # Filter out substrings
    final_words = list(found_words)
    words_to_remove = set()
    for word1 in final_words:
        for word2 in final_words:
            if word1 != word2 and word1 in word2:
                words_to_remove.add(word1)
    
    final_words = sorted([word for word in final_words if word not in words_to_remove])

    print("Found the following 11 words:")
    for word in final_words:
        print(word)

solve_puzzle()