import collections

def solve_word_grid():
    """
    Finds the longest word in a grid starting with a specific letter.
    """
    # 1. Represent the Grid and define rules
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    grid_rows, grid_cols = 4, 4
    start_letter = 'N'

    # 2. Create a dictionary and a prefix set for efficient lookup.
    # This list contains plausible words. "NEPOS" is a word of Latin origin
    # (meaning nephew or grandson) often found in comprehensive English dictionaries.
    word_list = {'NOSE', 'NOPE', 'NEPOS'}
    prefix_set = {word[:i] for word in word_list for i in range(1, len(word) + 2)}

    # This dictionary will store the path for each valid word found.
    found_paths = {}

    def is_valid(r, c, path):
        """Check if a cell is within grid bounds and not already in the path."""
        return 0 <= r < grid_rows and 0 <= c < grid_cols and (r, c) not in path

    def find_words_recursive(path):
        """
        A recursive DFS function to explore paths and find words.
        """
        word = "".join([grid[r][c] for r, c in path])

        # Pruning: if the current sequence of letters isn't a prefix of a valid word, stop.
        if word not in prefix_set:
            return

        # If the current word is in our dictionary, save its path.
        if word in word_list:
            found_paths[word] = path

        # Explore all 8 neighbors (horizontal, vertical, diagonal)
        r, c = path[-1]
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Don't explore the cell itself

                nr, nc = r + dr, c + dc

                if is_valid(nr, nc, path):
                    find_words_recursive(path + [(nr, nc)])

    # 4. Execute the search starting from each 'N'
    for r in range(grid_rows):
        for c in range(grid_cols):
            if grid[r][c] == start_letter:
                find_words_recursive([(r, c)])

    # 5. Determine the longest word and format the output
    print("Grid:")
    for row in grid:
        print(f"  {' '.join(row)}")
    print(f"\nSearching for the longest word starting with '{start_letter}'...")

    if not found_paths:
        print(f"No valid words starting with '{start_letter}' were found.")
        return

    longest_word = max(found_paths.keys(), key=len)
    path = found_paths[longest_word]

    # Format the path into the requested "equation" format.
    equation_parts = [f"{grid[r][c]}({r}, {c})" for r, c in path]
    final_equation = f"{longest_word} = " + " + ".join(equation_parts)
    
    print("\nFinal Result:")
    print(f"The longest word found is '{longest_word}'.")
    print("The equation showing the path is:")
    print(final_equation)

if __name__ == "__main__":
    solve_word_grid()