def find_word_path(grid, word):
    """
    Finds and prints the path of a given word in the grid.
    """
    if not word:
        print("Empty word provided.")
        return

    rows, cols = len(grid), len(grid[0])
    
    # This recursive helper function will perform the DFS
    def search(r, c, index, visited):
        # If we have found all characters of the word
        if index == len(word):
            return visited

        # Explore neighbors (including diagonals)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue

                nr, nc = r + dr, c + dc
                
                # Check if the neighbor is within bounds, not visited, and has the correct character
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited and grid[nr][nc] == word[index]:
                    visited.append((nr, nc))
                    result = search(nr, nc, index + 1, visited)
                    if result:
                        return result
                    # Backtrack if the path doesn't lead to a solution
                    visited.pop()
        return None

    # Start the search from any cell containing the first letter of the word
    start_char = word[0]
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == start_char:
                path_coords = search(r, c, 1, [(r, c)])
                if path_coords:
                    # Construct and print the final path string
                    path_letters = [grid[pr][pc] for pr, pc in path_coords]
                    print(" -> ".join(path_letters))
                    return
    
    print(f"The word '{word}' cannot be formed on the grid.")


# The 4x4 grid of letters
grid = [
    ['E', 'S', 'G', 'D'],
    ['K', 'P', 'O', 'J'],
    ['O', 'E', 'N', 'Y'],
    ['L', 'P', 'S', 'N']
]

# The target word to find
longest_word = "NOPES"

find_word_path(grid, longest_word)
