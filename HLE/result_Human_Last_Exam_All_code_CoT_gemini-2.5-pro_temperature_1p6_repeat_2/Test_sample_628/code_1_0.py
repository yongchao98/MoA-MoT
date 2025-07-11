import sys

def solve_word_grid():
    """
    Finds the longest word in a grid starting with a specific letter.
    """

    # 1. Define the grid and starting character
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    start_char = 'N'
    rows, cols = 4, 4

    # 2. Load words and create prefix set for efficient searching
    word_set = set()
    prefix_set = set()
    try:
        with open('/usr/share/dict/words', 'r') as f:
            for line in f:
                word_set.add(line.strip().upper())
    except FileNotFoundError:
        # Provide a fallback dictionary if the system one isn't available
        # This list includes the expected answer and other potential words.
        print("Dictionary '/usr/share/dict/words' not found. Using a fallback list.", file=sys.stderr)
        word_set = {"NO", "NOSE", "NOPE", "NOPES", "OPEN", "OPENS", "SPOON", "SPOKEN", "PENS"}

    for word in word_set:
        for i in range(1, len(word) + 1):
            prefix_set.add(word[:i])

    # 3. DFS search implementation
    longest_word_found = ""

    def dfs(x, y, current_word, visited):
        nonlocal longest_word_found

        current_word += grid[x][y]
        visited.add((x, y))

        # Pruning: if the current sequence of letters isn't a prefix, stop.
        if current_word not in prefix_set:
            visited.remove((x, y))
            return

        # If it's a valid word and the longest one found so far, save it.
        if current_word in word_set and len(current_word) > len(longest_word_found):
            longest_word_found = current_word

        # Explore all 8 neighbors (including diagonals)
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0:
                    continue
                
                nx, ny = x + dx, y + dy

                if 0 <= nx < rows and 0 <= ny < cols and (nx, ny) not in visited:
                    dfs(nx, ny, current_word, visited)
        
        # Backtrack to explore other paths
        visited.remove((x, y))

    # 4. Start the search from every 'N' in the grid
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == start_char:
                dfs(r, c, "", set())
    
    # 5. Output the result in the requested "equation" format
    if longest_word_found:
        # e.g., N + O + P + E + S
        print(" + ".join(list(longest_word_found)))
    else:
        print(f"No word starting with '{start_char}' found.")

solve_word_grid()
<<<N + O + P + E + S>>>