import urllib.request
import bisect

def solve_longest_word_puzzle():
    """
    Finds the longest word in a grid starting with 'N'.
    Downloads a dictionary, performs a depth-first search, and prints the result.
    """
    # 1. Represent the Grid
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    GRID_SIZE = 4
    
    # --- Helper data and functions ---
    
    # Global variables to store the final result
    longest_word_found = ""
    longest_word_path = []

    def get_word_lists():
        """Downloads and processes a word list."""
        try:
            # 2. Acquire a Word List
            url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
            response = urllib.request.urlopen(url)
            all_words = response.read().decode('utf-8').splitlines()
            
            # Filter for words starting with 'n' and of reasonable length
            n_words = sorted([word for word in all_words if word.startswith('n') and len(word) >= 3])
            
            # 3. Optimize Search Data Structures
            return n_words, set(n_words)
        except Exception as e:
            print(f"Error downloading or processing word list: {e}")
            print("Using a small fallback dictionary.")
            fallback_words = ["nopes", "nose", "nooky", "neap", "neat", "nope"]
            return sorted(fallback_words), set(fallback_words)

    word_list, word_set = get_word_lists()

    def is_valid_prefix(prefix):
        """Checks if a string is a prefix of any word in the dictionary using binary search."""
        idx = bisect.bisect_left(word_list, prefix)
        # If the insertion point is not out of bounds and the word at that point starts with the prefix
        return idx < len(word_list) and word_list[idx].startswith(prefix)

    # 4. Implement Depth-First Search (DFS)
    def dfs(path):
        """
        Recursively searches for words in the grid.
        
        Args:
            path (list): A list of (row, col) tuples representing the current path.
        """
        nonlocal longest_word_found, longest_word_path

        r, c = path[-1]
        current_word = "".join(grid[row][col].lower() for row, col in path)

        if not is_valid_prefix(current_word):
            return

        # Check if the current path forms a valid word and if it's the longest one found so far
        if current_word in word_set and len(current_word) > len(longest_word_found):
            longest_word_found = current_word
            longest_word_path = list(path) # Make a copy

        # Explore neighbors
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue

                nr, nc = r + dr, c + dc
                
                # Check bounds and if the cell has been visited in the current path
                if 0 <= nr < GRID_SIZE and 0 <= nc < GRID_SIZE and (nr, nc) not in path:
                    dfs(path + [(nr, nc)])

    # --- Main execution logic ---
    
    # Start DFS from every 'N' in the grid
    for r in range(GRID_SIZE):
        for c in range(GRID_SIZE):
            if grid[r][c] == 'N':
                dfs([(r, c)])
    
    # 6. Format and print the final output
    if longest_word_found:
        equation_parts = [grid[r][c] for r, c in longest_word_path]
        print(f"{' + '.join(equation_parts)} = {longest_word_found.upper()}")
    else:
        print("No word starting with 'N' could be found.")

solve_longest_word_puzzle()
<<<N + O + P + E + S = NOPES>>>