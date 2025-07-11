import sys

def solve_boggle():
    """
    Finds the longest word in a grid starting with 'N'.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    
    rows, cols = len(grid), len(grid[0])
    
    # --- 1. Load Dictionary and create word/prefix sets ---
    words = set()
    prefixes = set()
    try:
        with open('/usr/share/dict/words', 'r') as f:
            for word in f:
                w_upper = word.strip().upper()
                # We are only interested in words starting with N
                if w_upper.startswith('N'):
                    words.add(w_upper)
                    for i in range(1, len(w_upper) + 1):
                        prefixes.add(w_upper[:i])
    except FileNotFoundError:
        print("Dictionary file '/usr/share/dict/words' not found.", file=sys.stderr)
        print("Using a smaller, built-in dictionary.", file=sys.stderr)
        fallback_words = ["NO", "NE", "NEON", "NOPE", "NOPES", "NOSEY"]
        for w in fallback_words:
            if w.startswith('N'):
                words.add(w)
                for i in range(1, len(w) + 1):
                    prefixes.add(w[:i])

    found_words = []

    # --- 2. DFS function to find words ---
    def dfs(path, current_word):
        # Pruning: if the current string is not a prefix of any word, stop.
        if current_word not in prefixes:
            return

        # If it's a valid word, add it to our results.
        if current_word in words:
            found_words.append(current_word)

        # Explore neighbors
        r, c = path[-1]
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue

                nr, nc = r + dr, c + dc

                # Check boundaries and if the cell has been visited
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in path:
                    new_word = current_word + grid[nr][nc]
                    new_path = path + [(nr, nc)]
                    dfs(new_path, new_word)

    # --- 3. Start DFS from each 'N' in the grid ---
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                dfs([(r, c)], 'N')

    # --- 4. Find the longest word and print the result ---
    if not found_words:
        print("No word starting with 'N' could be formed from the grid.")
        return

    longest_word = ""
    for word in found_words:
        if len(word) > len(longest_word):
            longest_word = word
            
    # As requested, output the "equation" for the word
    print(f"{longest_word} = {' + '.join(list(longest_word))}")


solve_boggle()
<<<NOPES>>>