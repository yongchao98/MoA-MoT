import collections

def solve_word_grid():
    """
    Finds the longest word in a 4x4 grid starting with the letter 'N'.
    """
    # 1. Define the grid and word lists
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    
    # A curated word list. Adding a variety ensures the search is robust.
    word_list = [
        "NO", "NOES", "NOPE", "NOTE", "NOTES", "PENS", "POKE", "POKES", "SPOKE",
        "SPOKEN", "NSPOKE", "NESPOKE", "NESPOKES", "SEND", "SENT", "NEP"
    ]
    words = set(word_list)
    
    # 2. Create a prefix set for search optimization
    prefixes = set()
    for word in words:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    # 3. Implement the Depth-First Search (DFS)
    longest_word_info = {"word": "", "path": []}
    
    def is_valid_move(r, c, path_coords):
        """Checks if a move to (r, c) is within grid bounds and not already visited."""
        return 0 <= r < 4 and 0 <= c < 4 and (r, c) not in path_coords

    def find_longest_word_from(path):
        """Recursively search for words starting from a given path."""
        nonlocal longest_word_info
        
        # Form the current word from the path coordinates
        current_word = "".join([grid[r][c] for r, c in path])
        
        # Pruning: if the current sequence of letters isn't a valid prefix, stop.
        if current_word not in prefixes:
            return

        # If it's a valid word, check if it's the longest one found so far.
        if current_word in words:
            if len(current_word) > len(longest_word_info["word"]):
                longest_word_info["word"] = current_word
                longest_word_info["path"] = path

        # Continue the search from the current cell to its neighbors
        last_r, last_c = path[-1]
        
        # Explore all 8 directions (including diagonals)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Don't move to the same cell
                
                nr, nc = last_r + dr, last_c + dc
                
                if is_valid_move(nr, nc, path):
                    # Explore this new path
                    find_longest_word_from(path + [(nr, nc)])

    # 4. Start the search from every 'N' in the grid
    for r in range(4):
        for c in range(4):
            if grid[r][c] == 'N':
                # The initial path is just the starting cell
                find_longest_word_from([(r, c)])
    
    # 5. Print the result in the requested format
    final_word = longest_word_info["word"]
    final_path = longest_word_info["path"]
    
    if final_word:
        equation_parts = [grid[r][c] for r, c in final_path]
        equation_str = " + ".join(equation_parts)
        print(f"Found longest word: {final_word}")
        print(f"Path Coordinates: {final_path}")
        print(f"Equation: {equation_str} = {final_word}")
    else:
        print("No word starting with 'N' could be formed.")

if __name__ == '__main__':
    solve_word_grid()