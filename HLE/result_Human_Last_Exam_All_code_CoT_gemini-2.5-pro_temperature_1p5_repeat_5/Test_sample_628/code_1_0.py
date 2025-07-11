import sys
import urllib.request

def solve_word_grid():
    """
    Finds the longest word starting with 'N' in a 4x4 grid.
    The script downloads a word list, performs a depth-first search
    from all possible starting positions, and prints the longest word found.
    """
    # Step 1: Define the grid and starting character
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    start_char = 'N'
    rows, cols = 4, 4

    # Step 2: Download and prepare the dictionary
    try:
        # A common, comprehensive word list
        word_url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        print("Downloading word list...", file=sys.stderr)
        with urllib.request.urlopen(word_url) as response:
            word_data = response.read().decode('utf-8')
        words = set(word_data.split())
        print("Download complete.", file=sys.stderr)
    except Exception as e:
        print(f"Error: Could not download word list ({e}).", file=sys.stderr)
        print("Please check your internet connection and try again.", file=sys.stderr)
        # Using a small backup list for demonstration if download fails
        words = {"no", "nope", "nopes", "nose", "nosy"}
        print("Using a small backup dictionary.", file=sys.stderr)


    # Step 3: Create a prefix set for search optimization
    prefixes = {word[:i] for word in words for i in range(1, len(word) + 1)}

    # State variables for the search
    longest_word_found = ""
    best_path_coords = []

    def is_valid(r, c):
        """Check if coordinates are within the grid bounds."""
        return 0 <= r < rows and 0 <= c < cols

    def find_words_dfs(r, c, current_word, path_coords):
        """
        Recursively explores the grid to find words.
        Args:
            r (int): Current row.
            c (int): Current column.
            current_word (str): The word formed so far.
            path_coords (list): A list of (row, col) tuples visited in the current path.
        """
        nonlocal longest_word_found, best_path_coords

        # Append the current letter and position to the path
        current_word += grid[r][c]
        path_coords.append((r, c))

        # Pruning: If the current string is not a prefix of any word, stop this path.
        if current_word.lower() not in prefixes:
            return

        # Check if we have found a valid, longer word
        if current_word.lower() in words and len(current_word) > len(longest_word_found):
            longest_word_found = current_word
            best_path_coords = list(path_coords)

        # Explore all 8 neighbors (including diagonals)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the current cell itself
                
                nr, nc = r + dr, c + dc
                
                # If neighbor is valid and not already in our path, recurse
                if is_valid(nr, nc) and (nr, nc) not in path_coords:
                    find_words_dfs(nr, nc, current_word, list(path_coords))

    # Step 5: Start the DFS from each occurrence of the start_char
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == start_char:
                find_words_dfs(r, c, "", [])
    
    # Step 6: Print the final result
    if longest_word_found:
        print(f"The longest word found is '{longest_word_found}'.")
        # Construct the "equation" string from the saved path
        path_letters = [grid[r][c] for r, c in best_path_coords]
        print("It is formed by the following sequence of letters:")
        print(" -> ".join(path_letters))
    else:
        print(f"No word starting with '{start_char}' could be found in the grid.")

if __name__ == "__main__":
    solve_word_grid()
