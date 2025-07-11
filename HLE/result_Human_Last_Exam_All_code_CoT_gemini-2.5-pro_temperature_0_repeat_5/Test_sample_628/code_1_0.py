import collections

def solve_boggle():
    """
    Finds the longest word starting with 'N' in a 4x4 grid.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # A small dictionary for this specific problem.
    # In a real-world scenario, this would be loaded from a large word file.
    word_list = {"NO", "NOPE", "NOPES", "NOSE", "SPOKEN"}
    
    # Create a prefix set for efficient pruning during the search.
    prefixes = set()
    for word in word_list:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    longest_word = ""

    def is_valid(x, y):
        return 0 <= x < 4 and 0 <= y < 4

    def find_words_recursive(x, y, current_word, visited):
        nonlocal longest_word

        # Append current letter and mark as visited
        current_word += grid[x][y]
        visited.add((x, y))

        # Prune the search if the current string is not a valid prefix
        if current_word not in prefixes:
            # Backtrack
            visited.remove((x, y))
            return

        # If it's a valid word and longer than the current longest, update it
        if current_word in word_list and len(current_word) > len(longest_word):
            longest_word = current_word

        # Explore all 8 neighbors
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0:
                    continue
                
                nx, ny = x + dx, y + dy

                if is_valid(nx, ny) and (nx, ny) not in visited:
                    find_words_recursive(nx, ny, current_word, visited)
        
        # Backtrack: un-mark the cell for other paths
        visited.remove((x, y))

    # Find all starting positions for 'N'
    start_positions = []
    for r in range(4):
        for c in range(4):
            if grid[r][c] == 'N':
                start_positions.append((r, c))

    # Start the search from each 'N'
    for r, c in start_positions:
        find_words_recursive(r, c, "", set())

    print(longest_word)

solve_boggle()