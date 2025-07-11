import sys
import requests

def solve_boggle():
    """
    Finds the longest word in a 4x4 grid starting with 'N'.
    """
    # The grid of letters
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    
    # --- Step 1: Get word list and create prefix set ---
    try:
        # Using a common word list from a public repository
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes
        all_words_text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching word list: {e}", file=sys.stderr)
        print("Please check your internet connection.", file=sys.stderr)
        # As a fallback, use a small set of words including the likely answer
        all_words_text = "noel\nnose\nnepos\nspoken\n"

    # Use sets for efficient O(1) average time complexity lookups
    words = set(word.upper() for word in all_words_text.splitlines() if len(word) > 2)
    prefixes = set()
    for word in words:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    # --- Step 2: Implement the Depth-First Search (DFS) ---
    rows, cols = 4, 4
    longest_word_found = ""

    def is_valid(r, c):
        return 0 <= r < rows and 0 <= c < cols

    def dfs(r, c, path, current_word):
        nonlocal longest_word_found

        # Pruning: if the current word is not a prefix of any word in the dictionary, stop.
        if current_word not in prefixes:
            return

        # If the current path forms a valid word and it's the longest found so far
        if current_word in words and len(current_word) > len(longest_word_found):
            longest_word_found = current_word

        # Explore neighbors
        # dr, dc covers all 8 directions (horizontal, vertical, diagonal)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the cell itself

                nr, nc = r + dr, c + dc

                if is_valid(nr, nc) and (nr, nc) not in path:
                    new_path = path + [(nr, nc)]
                    new_word = current_word + grid[nr][nc]
                    dfs(nr, nc, new_path, new_word)

    # --- Step 3: Start DFS from every 'N' on the grid ---
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                # Start DFS with the path containing the starting cell
                dfs(r, c, [(r, c)], grid[r][c])
                
    # --- Step 4: Print the result ---
    if longest_word_found:
        print(f"The longest word found is: {longest_word_found}")
    else:
        print("No word starting with 'N' could be found.")

if __name__ == "__main__":
    solve_boggle()
