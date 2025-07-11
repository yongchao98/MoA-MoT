import requests
import io

def solve_word_search():
    """
    Finds the 11 hidden words in the provided grid based on the given constraints.
    """
    # Step 1: Define the grid
    grid = [
        "DESERTGFSG",
        "EEHAJWNLPS",
        "ILONSHIAED",
        "FKUWZEKMEU",
        "ICLHONNECO",
        "RIDKQEISHL",
        "TFIDMVHCLC",
        "EYADRETSEY",
        "PERLBRFPST",
        "BREATHLESS",
    ]
    grid_lower = [row.lower() for row in grid]
    rows, cols = len(grid_lower), len(grid_lower[0])

    # Step 2: Fetch a comprehensive English word list
    try:
        # Using a large, standard word list for better coverage
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        all_words = {line.strip().lower() for line in response.text.splitlines()}
    except requests.exceptions.RequestException as e:
        print(f"Could not download word list: {e}. The script cannot continue.")
        return

    # Filter for words of at least 6 letters for efficiency
    words_to_check = {word for word in all_words if len(word) >= 6}

    # Step 3: Search the grid in all 8 directions
    # Directions: (dr, dc) for N, NE, E, SE, S, SW, W, NW
    directions = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]
    found_words_with_info = set()

    for r_start in range(rows):
        for c_start in range(cols):
            for dr, dc in directions:
                word = ""
                r, c = r_start, c_start
                while 0 <= r < rows and 0 <= c < cols:
                    word += grid_lower[r][c]
                    if len(word) >= 6 and word in words_to_check:
                        # Store the found word along with its starting position
                        found_words_with_info.add((word, r_start, c_start))
                    r += dr
                    c += dc

    # Step 4: Filter out words that are substrings of other found words
    all_found_strings = {info[0] for info in found_words_with_info}
    final_words_info = set()

    for word_info in found_words_with_info:
        word_str = word_info[0]
        is_substring = False
        for other_word_str in all_found_strings:
            if word_str != other_word_str and word_str in other_word_str:
                is_substring = True
                break
        if not is_substring:
            final_words_info.add(word_info)

    # Step 5: Sort the final list by starting position (row, then column)
    sorted_words_info = sorted(list(final_words_info), key=lambda item: (item[1], item[2]))
    
    # Step 6: Print the final ordered list of words
    final_words = [info[0] for info in sorted_words_info]
    
    print("The 11 found words in order are:")
    for word in final_words:
        print(word)

if __name__ == '__main__':
    solve_word_search()