import requests

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
    
    # 1. Get a word list and prepare sets for fast lookups.
    # The dictionary file has lowercase words, so we convert them to uppercase to match the grid.
    try:
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        words = {word.upper() for word in response.text.split()}
    except requests.exceptions.RequestException:
        print("Could not download word list, using a small fallback list.")
        words = {"NO", "NOSE", "SPOKEN", "SPOON", "PENS", "GOSPEL", "NOSEY", "SNOOPY"}

    # Create a set of all possible prefixes for search pruning.
    prefixes = {word[:i] for word in words for i in range(1, len(word) + 1)}

    # 2. Initialize variables to store the result.
    rows, cols = len(grid), len(grid[0])
    longest_word_found = ""

    # 3. Implement the recursive Depth-First Search (DFS).
    def search(r, c, current_word, path):
        nonlocal longest_word_found

        # Form the new word and path
        letter = grid[r][c]
        new_word = current_word + letter
        new_path = path + [(r, c)]

        # Prune the search if the current sequence of letters is not a valid prefix.
        if new_word not in prefixes:
            return

        # If the new word is a valid word in our dictionary, check if it's the longest one found so far.
        if new_word in words:
            if len(new_word) > len(longest_word_found):
                longest_word_found = new_word

        # Explore all 8 neighbors.
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the current cell itself.

                nr, nc = r + dr, c + dc

                # Check if the neighbor is within grid boundaries and has not been visited yet in the current path.
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in new_path:
                    search(nr, nc, new_word, new_path)

    # 4. Start the search from every 'N' on the grid.
    start_char = 'N'
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == start_char:
                search(r, c, "", [])

    # 5. Print the final answer.
    if longest_word_found:
        print(f"The longest word starting with 'N' is: {longest_word_found}")
        # The prompt mentioned "final equation", which is likely a template error.
        # I will output the letters of the word as requested.
        print("The letters in the word are:")
        for char in longest_word_found:
            print(char)
    else:
        print(f"No word starting with '{start_char}' could be found.")

if __name__ == '__main__':
    solve_boggle()
