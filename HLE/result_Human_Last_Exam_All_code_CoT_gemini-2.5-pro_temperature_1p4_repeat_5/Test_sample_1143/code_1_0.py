import urllib.request
import sys

def solve_word_search_and_find_album():
    """
    Solves the word search puzzle to find 11 hidden words, then explains
    how these words lead to the final answer.
    """
    # Step 1: Download a list of English words.
    try:
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        print("Downloading English dictionary...", file=sys.stderr)
        with urllib.request.urlopen(url, timeout=15) as response:
            all_words = {line.decode('utf-8').strip().lower() for line in response}
        print("Dictionary downloaded.", file=sys.stderr)
        word_list = {word for word in all_words if len(word) >= 6}
    except Exception as e:
        print(f"Error: Could not download the dictionary from {url}.", file=sys.stderr)
        print(f"Please check your internet connection. Reason: {e}", file=sys.stderr)
        return

    # Step 2: Define the puzzle grid.
    grid_str = """
    DESERTGFSG
    EEHAJWNLPS
    ILONSHIAED
    FKUWZEKMEU
    ICLHONNECO
    RIDKQEISHL
    TFIDMVHCLC
    EYADRETSEY
    PERLBRFPST
    BREATHLESS
    """
    grid = [list(line.lower()) for line in grid_str.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])

    # Step 3: Search the grid for words.
    found_words = set()
    directions = [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]

    for r in range(rows):
        for c in range(cols):
            for dr, dc in directions:
                word = ""
                cr, cc = r, c
                while 0 <= cr < rows and 0 <= cc < cols:
                    word += grid[cr][cc]
                    if len(word) >= 6 and word in word_list:
                        found_words.add(word)
                    cr += dr
                    cc += dc
    
    # Step 4: Filter out words that are substrings of other found words.
    # This ensures that if "rain" and "rainbow" are found, only "rainbow" is kept.
    final_words = []
    found_list_sorted = sorted(list(found_words), key=len, reverse=True)

    for word in found_list_sorted:
        is_substring = False
        # Check if the current word is a substring of any word already in our final list.
        # Since we sorted by length descending, we only need to check against words
        # already added to final_words.
        for other_word in final_words:
            if word in other_word:
                is_substring = True
                break
        if not is_substring:
            final_words.append(word)
    
    # Step 5: Sort the final list of 11 words alphabetically and print them.
    final_words.sort()
    
    print("The 11 found words, sorted alphabetically, are:")
    for word in final_words:
        print(word)

if __name__ == '__main__':
    solve_word_search_and_find_album()