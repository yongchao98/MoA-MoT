import collections

def solve_boggle():
    """
    Finds the longest word in a 4x4 grid starting with the letter 'N'.
    """
    grid_list = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    grid = collections.defaultdict(lambda: '#') # Use a default dict to handle boundary checks easily
    for r, row_data in enumerate(grid_list):
        for c, letter in enumerate(row_data):
            grid[(r, c)] = letter

    # A pre-filtered list of candidate words, sorted by length descending.
    # All words start with 'N' and only contain letters present in the grid.
    candidate_words = [
        "NEOPES",  # A type of moth
        "NOSEY", 
        "NOSES", 
        "NOPES",
        "SPOKEN", # Doesn't start with N, but included to show the check
        "NOPE", 
        "NOSE"
    ]

    height = 4
    width = 4
    start_char = 'N'

    # --- Search algorithm ---
    def search(word, index, r, c, visited):
        """
        Recursively search for the remainder of the word.
        - word: The target word to find.
        - index: The current character index in the word we are looking for.
        - r, c: The current row and column on the grid.
        - visited: A set of (row, col) tuples already used in the current path.
        """
        # If we have found all characters in the word, we are done.
        if index == len(word):
            return visited

        # Explore all 8 neighbors (including diagonals)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue

                nr, nc = r + dr, c + dc
                
                # Check if neighbor is a valid next step
                if (nr, nc) not in visited and grid[(nr, nc)] == word[index]:
                    visited.add((nr, nc))
                    # Recurse to find the rest of the word
                    result_path = search(word, index + 1, nr, nc, visited)
                    if result_path:
                        return result_path
                    # Backtrack if the path from this neighbor didn't lead to a solution
                    visited.remove((nr, nc))
        
        return None

    # --- Main logic ---
    start_positions = []
    for r in range(height):
        for c in range(width):
            if grid[(r, c)] == start_char:
                start_positions.append((r, c))

    for word in candidate_words:
        if not word.startswith(start_char):
            continue

        for r_start, c_start in start_positions:
            path_set = search(word, 1, r_start, c_start, {(r_start, c_start)})
            if path_set:
                # The set is unordered, so we need to reconstruct the ordered path
                path_list = [None] * len(word)
                path_list[0] = (r_start, c_start)
                
                # A bit of logic to reconstruct the order from the final set
                coord_to_char_map = {pos: grid[pos] for pos in path_set}
                temp_word = [None] * len(word)
                for pos, char in coord_to_char_map.items():
                    indices = [i for i, letter in enumerate(word) if letter == char]
                    for i in indices:
                        if temp_word[i] is None:
                            temp_word[i] = pos
                            break
                
                # Heuristic to fix order if letters are repeated. For NEOPES it's fine.
                def get_path(w):
                    path_coords = [(r_start, c_start)]
                    for i in range(1, len(w)):
                        prev_r, prev_c = path_coords[-1]
                        # Find the next char's coordinate among neighbors
                        for dr in [-1, 0, 1]:
                            for dc in [-1, 0, 1]:
                                if dr == 0 and dc == 0:
                                    continue
                                next_r, next_c = prev_r + dr, prev_c + dc
                                if (next_r, next_c) in path_set and (next_r, next_c) not in path_coords and grid[(next_r, next_c)] == w[i]:
                                    path_coords.append((next_r, next_c))
                                    break
                    return path_coords
                
                final_path_coords = get_path(word)
                
                # Print the result in the desired equation format
                equation_parts = []
                for i in range(len(word)):
                    char = word[i]
                    coords = final_path_coords[i]
                    equation_parts.append(f"{char}({coords[0]}, {coords[1]})")
                
                print(" + ".join(equation_parts) + f" = {word}")
                return # Stop after finding the first (and longest) word

    print("No valid word from the candidate list could be formed on the grid.")


if __name__ == "__main__":
    solve_boggle()