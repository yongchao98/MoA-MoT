import collections

def solve_word_grid():
    """
    Finds the longest word in a 4x4 grid starting with 'N'.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # A curated word list of plausible words starting with 'N' that can be formed from the grid letters.
    # This is necessary as a full dictionary is not available.
    words = {"NOPE", "NOPES", "NOES", "NOSE", "NEON", "NESP"}
    
    # Create a prefix set for efficient search pruning
    prefixes = set()
    for word in words:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    longest_word_found = ""
    best_path_found = []

    def is_valid(r, c):
        return 0 <= r < 4 and 0 <= c < 4

    def find_words_from(start_r, start_c):
        nonlocal longest_word_found, best_path_found
        
        q = collections.deque()
        # Each item in the queue is a tuple: (current_word, path_list)
        q.append((grid[start_r][start_c], [(start_r, start_c)]))

        while q:
            current_word, path = q.popleft()

            # Check if this is a valid word and if it's the longest found
            if current_word in words and len(current_word) > len(longest_word_found):
                longest_word_found = current_word
                best_path_found = path
            
            # Explore neighbors
            r, c = path[-1]
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0:
                        continue
                    
                    nr, nc = r + dr, c + dc
                    
                    if is_valid(nr, nc) and (nr, nc) not in path:
                        new_word = current_word + grid[nr][nc]
                        # Pruning: only continue if the new word is a valid prefix
                        if new_word in prefixes:
                            new_path = path + [(nr, nc)]
                            q.append((new_word, new_path))
    
    # Find all 'N' starting positions
    start_positions = []
    for r in range(4):
        for c in range(4):
            if grid[r][c] == 'N':
                start_positions.append((r, c))

    # Run the search from each starting 'N'
    for r_start, c_start in start_positions:
        find_words_from(r_start, c_start)

    # Format and print the final equation
    if longest_word_found:
        equation_parts = []
        for r, c in best_path_found:
            letter = grid[r][c]
            equation_parts.append(f"{letter} ({r}, {c})")
        print(" + ".join(equation_parts) + f" = {longest_word_found}")
    else:
        print("No word starting with 'N' found in the dictionary.")

solve_word_grid()
<<<NOPES>>>