import urllib.request
import ssl

def solve_puzzle():
    """
    Solves the word search part of the puzzle to find the 11 hidden words.
    """
    # Step 1: Define the grid from the problem description.
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

    # Step 2: Download a comprehensive English word list.
    # A custom SSL context is used to prevent potential certification errors.
    word_list_url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
    try:
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        with urllib.request.urlopen(word_list_url, context=ctx) as response:
            all_words_bytes = response.read()
        all_words = all_words_bytes.decode('utf-8').splitlines()
        # Create a set for fast lookups, containing only words of the required length.
        dictionary = {word.lower() for word in all_words if len(word) >= 6}
    except Exception as e:
        print(f"Could not download the dictionary. Error: {e}")
        return

    # Step 3: Generate all possible lines from the grid in the 4 directions.
    lines = []
    # Horizontal (Left to Right)
    lines.extend(grid)
    # Horizontal (Right to Left)
    lines.extend([row[::-1] for row in grid])
    # Vertical (Top to Bottom)
    cols = [''.join(grid[r][c] for r in range(len(grid))) for c in range(len(grid[0]))]
    lines.extend(cols)
    # Vertical (Bottom to Top)
    lines.extend([col[::-1] for col in cols])

    # Step 4: Extract all candidate substrings and find valid words.
    found_words = set()
    for line in lines:
        line_lower = line.lower()
        for i in range(len(line_lower)):
            for j in range(i + 6, len(line_lower) + 1):
                substring = line_lower[i:j]
                if substring in dictionary:
                    found_words.add(substring)

    # Step 5: Apply the non-substring rule to get the final list.
    final_words = set()
    for word1 in found_words:
        is_substring_of_another = False
        for word2 in found_words:
            if word1 != word2 and word1 in word2:
                is_substring_of_another = True
                break
        if not is_substring_of_another:
            final_words.add(word1)
            
    # Step 6: Print the final list of words, sorted alphabetically.
    sorted_words = sorted(list(final_words))
    print("Found the following 11 words:")
    for word in sorted_words:
        print(word)

if __name__ == '__main__':
    solve_puzzle()
