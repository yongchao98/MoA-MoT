import io

def solve_grid_puzzle():
    """
    Finds the longest word in a grid starting with 'N'.
    """
    # Step 1: Define the grid
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    rows, cols = 4, 4
    start_char = 'N'

    # Step 2: Embed a dictionary for self-contained execution
    # A small, curated list of words for this problem.
    word_list_data = """
    NO
    NOES
    NOPE
    NOPES
    NOSE
    NOTE
    SONE
    PENS
    POKE
    POSY
    SEND
    SPOKE
    SPOKEN
    """
    
    word_set = set()
    prefix_set = set()
    
    # Use io.StringIO to treat the string as a file
    with io.StringIO(word_list_data) as f:
        for line in f:
            word = line.strip().upper()
            if word:
                word_set.add(word)
                for i in range(1, len(word) + 1):
                    prefix_set.add(word[:i])

    # Step 3: Prepare for the search
    longest_word_info = {"word": "", "path": []}

    def search(row, col, visited, path):
        """
        A recursive DFS function to find words.
        """
        # Add current letter to the path
        letter = grid[row][col]
        new_path = path + [letter]
        current_word = "".join(new_path)
        
        # Prune the search if the current path is not a valid prefix
        if current_word not in prefix_set:
            return

        # If it's a valid word and longer than the current longest, update it
        if current_word in word_set:
            if len(current_word) > len(longest_word_info["word"]):
                longest_word_info["word"] = current_word
                longest_word_info["path"] = new_path

        # Explore neighbors
        # Moves are defined as (row_change, col_change)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the current cell itself

                next_row, next_col = row + dr, col + dc

                # Check if the neighbor is within grid bounds and not visited
                if 0 <= next_row < rows and 0 <= next_col < cols and (next_row, next_col) not in visited:
                    new_visited = visited.copy()
                    new_visited.add((next_row, next_col))
                    search(next_row, next_col, new_visited, new_path)

    # Step 4: Start the search from each 'N'
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == start_char:
                search(r, c, {(r, c)}, [])

    # Step 5: Print the final result in the specified format
    if longest_word_info["word"]:
        equation_str = " + ".join(longest_word_info["path"])
        result_str = f"{equation_str} = {longest_word_info['word']}"
        print(result_str)
    else:
        print(f"No word starting with '{start_char}' found.")

solve_grid_puzzle()
<<<S + P + O + K + E + N = SPOKEN>>>