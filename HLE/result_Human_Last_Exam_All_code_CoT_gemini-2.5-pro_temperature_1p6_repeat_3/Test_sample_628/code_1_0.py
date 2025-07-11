import collections

def solve_boggle():
    """
    Finds the longest word in a 4x4 grid starting with 'N'.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # Using a curated dictionary for a self-contained script.
    # A larger dictionary like /usr/share/dict/words would be more comprehensive.
    dictionary = {
        "NO", "NOPE", "NOPES", "NOES", "GOES", "PENS", "SEND",
        "SON", "SONY", "POSE", "POKE", "POKES"
    }
    
    # Create a prefix set for efficient pruning during the search
    prefix_set = set()
    for word in dictionary:
        for i in range(1, len(word) + 1):
            prefix_set.add(word[:i])

    rows, cols = 4, 4
    longest_word_info = {"word": "", "path": []}

    def dfs(x, y, current_path, current_word):
        """
        Recursively search for words in the grid.
        
        Args:
            x (int): current row index.
            y (int): current column index.
            current_path (list): list of (row, col) tuples in the current path.
            current_word (str): the word formed by the current path.
        """
        nonlocal longest_word_info
        
        # Pruning step: if the current word is not a valid prefix, stop this path.
        if current_word not in prefix_set:
            return

        # Check if we have found a valid, longer word
        if current_word in dictionary:
            if len(current_word) > len(longest_word_info["word"]):
                longest_word_info["word"] = current_word
                longest_word_info["path"] = current_path
        
        # Explore neighbors
        for i in range(-1, 2):
            for j in range(-1, 2):
                if i == 0 and j == 0:
                    continue
                
                nx, ny = x + i, y + j

                # Check if neighbor is within bounds and not already in the path
                if 0 <= nx < rows and 0 <= ny < cols and (nx, ny) not in current_path:
                    dfs(nx, ny, current_path + [(nx, ny)], current_word + grid[nx][ny])

    # Start the DFS from every 'N' in the grid
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                dfs(r, c, [(r, c)], 'N')
    
    # Format and print the final result
    if longest_word_info["word"]:
        path_letters = [grid[r][c] for r, c in longest_word_info["path"]]
        equation = " + ".join(path_letters)
        print(f"{equation} = {longest_word_info['word']}")
    else:
        print("No word starting with 'N' was found.")

solve_boggle()