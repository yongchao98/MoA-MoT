import urllib.request
import sys

def solve_boggle():
    """
    Solves the Boggle grid puzzle to find the longest word starting with a specific letter.
    This function is self-contained and will download its own dictionary.
    """
    # The 4x4 grid of letters for the puzzle.
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    # The required starting letter for the word.
    start_char = 'N'
    rows, cols = len(grid), len(grid[0])

    # --- Dictionary Loading ---
    print("Loading dictionary... This might take a moment.")
    words = set()
    prefixes = set()
    # A reliable URL for a large English word list.
    url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
    try:
        # Open the URL and read the content.
        with urllib.request.urlopen(url) as response:
            # Decode the text and split into words.
            text = response.read().decode('utf-8')
            # Convert all words to uppercase to match the grid.
            word_list = [word.upper() for word in text.splitlines()]
            
            # Populate the sets for fast lookups.
            words = set(word_list)
            for word in word_list:
                for i in range(1, len(word) + 1):
                    prefixes.add(word[:i])
    except Exception as e:
        print(f"Error: Could not download or process the dictionary from {url}", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
        print("Please check your internet connection and try again.", file=sys.stderr)
        sys.exit(1)
    print("Dictionary loaded successfully.")

    # --- Search Algorithm (DFS) ---
    # Using a list to hold the result so it's mutable across recursive calls.
    longest_word_found = [""]

    def is_valid(r, c, path):
        """Checks if a cell is within grid bounds and not already in the current path."""
        return 0 <= r < rows and 0 <= c < cols and (r, c) not in path

    def dfs(r, c, path, current_word):
        """
        Performs a depth-first search from cell (r, c) to find words.
        
        Args:
            r (int): Current row index.
            c (int): Current column index.
            path (set): A set of (row, col) tuples of visited cells for the current word.
            current_word (str): The word formed so far along the current path.
        """
        # Pruning: If the current string is not a prefix of any valid word, stop this search path.
        if current_word not in prefixes:
            return

        # Check if the current path forms a valid word and if it's the longest one found so far.
        if current_word in words and len(current_word) > len(longest_word_found[0]):
            longest_word_found[0] = current_word

        # Explore all 8 neighbors (horizontal, vertical, and diagonal).
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                # Skip the current cell itself.
                if dr == 0 and dc == 0:
                    continue

                next_r, next_c = r + dr, c + dc

                # If the neighbor is a valid next step, continue the search from there.
                if is_valid(next_r, next_c, path):
                    new_word = current_word + grid[next_r][next_c]
                    dfs(next_r, next_c, path | {(next_r, next_c)}, new_word)

    # --- Main Execution ---
    # Find all starting positions for the given start_char.
    start_positions = []
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == start_char:
                start_positions.append((r, c))

    print(f"\nSearching for the longest word starting with '{start_char}'...")
    # Start the DFS from each identified starting position.
    for r_start, c_start in start_positions:
        dfs(r_start, c_start, {(r_start, c_start)}, start_char)

    # --- Final Output ---
    final_word = longest_word_found[0]
    if final_word:
        print("\n--- Result ---")
        print(f"The longest word that can be formed is: {final_word}")
    else:
        print(f"\nNo word starting with '{start_char}' could be found in the grid.")

if __name__ == "__main__":
    solve_boggle()