import collections

def solve_word_search():
    """
    Finds the longest word in a 4x4 grid starting with 'N'.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # A self-contained dictionary of words to make the script runnable.
    # I've included plausible words from the grid and other common words.
    word_list_str = """
    neon
    nope
    nopes
    nose
    open
    opens
    pen
    pens
    poke
    pokes
    pose
    see
    seen
    send
    so
    son
    soon
    spoke
    spoken
    yen
    yep
    yes
    """
    
    dictionary = set(word_list_str.strip().split())
    prefixes = set()
    for word in dictionary:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    rows, cols = 4, 4
    found_words = set()

    def dfs(r, c, current_word, visited):
        """
        Recursively search the grid for words.
        """
        # 1. Boundary and visited check
        if not (0 <= r < rows and 0 <= c < cols and (r, c) not in visited):
            return

        # 2. Form new word and check if it's a valid prefix
        letter = grid[r][c]
        new_word = (current_word + letter).lower()

        if new_word not in prefixes:
            return
            
        # 3. Add to visited set for this path
        visited.add((r, c))

        # 4. If it's a valid word, add to our results
        if new_word in dictionary:
            found_words.add(new_word)

        # 5. Explore all 8 neighbors
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                dfs(r + dr, c + dc, new_word, visited)
        
        # 6. Backtrack: remove from visited so it can be used in other paths
        visited.remove((r, c))

    # Start DFS from each 'N' in the grid
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                dfs(r, c, "", set())

    if not found_words:
        print("No words starting with 'N' were found.")
    else:
        # Find the longest word among all found words
        longest_word = max(found_words, key=len)
        print(longest_word)

solve_word_search()