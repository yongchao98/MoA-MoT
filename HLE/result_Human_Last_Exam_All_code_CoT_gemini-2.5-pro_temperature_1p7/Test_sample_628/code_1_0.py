import urllib.request
import sys

def solve_word_grid():
    """
    Finds the longest word in a grid starting with a specific letter.
    """
    # 1. Represent the Grid
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # 2. Get a Dictionary
    try:
        print("Downloading dictionary file...", file=sys.stderr)
        # A common, comprehensive word list
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        with urllib.request.urlopen(url) as response:
            # Decode, split into lines, and convert to an uppercase set for fast lookups
            text = response.read().decode('utf-8')
            dictionary = {word.strip().upper() for word in text.splitlines()}
        print("Dictionary loaded.", file=sys.stderr)
    except Exception as e:
        print(f"Error: Could not download or process dictionary file: {e}", file=sys.stderr)
        print("Please check your internet connection and try again.", file=sys.stderr)
        return

    # 4. Explore Paths (Depth-First Search)
    rows, cols = len(grid), len(grid[0])
    longest_word = ""

    # This recursive function will explore all paths from a given cell
    def dfs(r, c, current_word, visited_coords):
        nonlocal longest_word

        # Add current state to the path and word
        current_word += grid[r][c]
        visited_coords.add((r, c))

        # Check if the current word is valid and the new longest
        if current_word in dictionary and len(current_word) > len(longest_word):
            longest_word = current_word

        # Explore all 8 neighbors (horizontal, vertical, diagonal)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the cell itself

                nr, nc = r + dr, c + dc

                # Check if the neighbor is within bounds and not already visited
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited_coords:
                    # Note: We pass a copy of visited_coords to the recursive call
                    # so that different search branches don't affect each other.
                    dfs(nr, nc, current_word, visited_coords.copy())

    # 3. Find Starting Points and Initiate Search
    start_letter = 'N'
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == start_letter:
                # Start a new search from this cell
                dfs(r, c, "", set())

    # 6. Format the Output
    if longest_word:
        # Print each letter of the word separated by " + "
        print(" + ".join(list(longest_word)))
    else:
        print(f"No valid word starting with '{start_letter}' was found.")

if __name__ == "__main__":
    solve_word_grid()