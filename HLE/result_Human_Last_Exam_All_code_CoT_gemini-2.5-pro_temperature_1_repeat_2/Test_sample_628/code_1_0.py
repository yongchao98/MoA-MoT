import urllib.request

def solve_boggle():
    """
    Finds the longest word in a grid starting with a specific letter.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    start_char = 'N'
    
    # --- Dictionary and Prefix Set Setup ---
    try:
        # Fetch a word list from a public repository
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        with urllib.request.urlopen(url) as response:
            word_list = response.read().decode('utf-8').splitlines()
        print("Successfully loaded dictionary from the web.")
    except Exception as e:
        print(f"Could not load dictionary from the web ({e}). Using a small backup list.")
        # Backup list in case of network failure
        word_list = ["no", "nope", "nopes", "nose", "spoken", "poison", "poisons", "nyse", "neop"]

    # Filter words and create efficient lookup sets for words and prefixes
    words = {word.upper() for word in word_list if len(word) >= 3}
    prefixes = set()
    for word in words:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])
            
    # --- DFS Solver ---
    rows, cols = len(grid), len(grid[0])
    longest_word = ""

    def solve_recursive(r, c, current_word, visited):
        nonlocal longest_word
        
        # Form the new word by adding the current letter
        current_word += grid[r][c]
        
        # 1. Pruning: If the current combination of letters isn't a prefix of any word, stop.
        if current_word not in prefixes:
            return

        # 2. Candidate Found: If it's a valid word and longer than our current best, update it.
        if current_word in words and len(current_word) > len(longest_word):
            longest_word = current_word

        # 3. Exploration: Recurse on all valid neighbors
        # Directions: Up, Down, Left, Right, and all 4 diagonals
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue # Skip the cell itself

                nr, nc = r + dr, c + dc

                # Check boundaries and if the cell has been visited
                if 0 <= nr < rows and 0 <= nc < cols and not visited[nr][nc]:
                    visited[nr][nc] = True
                    solve_recursive(nr, nc, current_word, visited)
                    visited[nr][nc] = False # Backtrack

    # --- Main Execution ---
    # Iterate over the grid to find all possible starting points
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == start_char:
                # For each start, create a fresh 'visited' grid
                visited = [[False for _ in range(cols)] for _ in range(rows)]
                visited[r][c] = True
                solve_recursive(r, c, "", visited)

    print(f"The longest word found starting with '{start_char}' is:")
    print(longest_word)

solve_boggle()
<<<NOPES>>>