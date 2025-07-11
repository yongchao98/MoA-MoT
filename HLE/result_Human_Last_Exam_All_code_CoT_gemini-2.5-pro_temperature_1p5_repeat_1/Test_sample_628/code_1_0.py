import urllib.request

def solve_word_grid():
    """
    Finds the longest word in a grid starting with a specific letter.
    """
    # The grid of letters
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    
    # Global variable to store the result
    longest_word_found = ""

    def get_dictionary():
        """
        Fetches a dictionary from a URL and creates word and prefix sets.
        Includes a fallback dictionary in case of network errors.
        """
        try:
            url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
            print("Fetching dictionary from the web...")
            with urllib.request.urlopen(url, timeout=5) as response:
                text = response.read().decode('utf-8')
            print("Dictionary loaded.")
            words = text.splitlines()
            word_set = {word.upper() for word in words}
        except Exception as e:
            print(f"Could not fetch dictionary from web ({e}). Using a smaller fallback list.")
            fallback_words = {"NOSE", "NOSEY", "NEAT", "NOTE", "SONY", "SPOKE", "SPOKEN", "LENS", "PENS"}
            word_set = {w.upper() for w in fallback_words}
            
        prefix_set = {word[:i] for word in word_set for i in range(1, len(word) + 1)}
        return word_set, prefix_set

    word_set, prefix_set = get_dictionary()
    rows, cols = len(grid), len(grid[0])

    def dfs(path):
        """
        Recursive Depth-First Search function to find words.
        'path' is a list of (row, col) tuples.
        """
        nonlocal longest_word_found
        
        # Construct the current word from the path of coordinates
        current_word = "".join(grid[r][c] for r, c in path)
        
        # --- Optimization Step ---
        # If the current string is not a prefix of any valid word, stop this path.
        if current_word not in prefix_set:
            return

        # --- Found a valid word ---
        # If the word is in our dictionary and is longer than the current longest, update it.
        if current_word in word_set and len(current_word) > len(longest_word_found):
            longest_word_found = current_word

        # --- Recursive Step ---
        # Explore all 8 neighbors of the last cell in the path.
        r, c = path[-1]
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the cell itself
                
                nr, nc = r + dr, c + dc
                
                # Check if the neighbor is within grid bounds and not already in the path.
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in path:
                    # Continue the search from the neighbor
                    dfs(path + [(nr, nc)])

    # --- Main Loop ---
    # Start the DFS from every cell containing 'N'.
    print("Searching for the longest word starting with 'N'...")
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                dfs([(r, c)])
    
    # --- Final Output ---
    if longest_word_found:
        print(f"\nThe longest word found is: {longest_word_found}")
    else:
        print("\nNo valid word starting with 'N' was found.")

# Run the solver
solve_word_grid()
<<<NOSEY>>>