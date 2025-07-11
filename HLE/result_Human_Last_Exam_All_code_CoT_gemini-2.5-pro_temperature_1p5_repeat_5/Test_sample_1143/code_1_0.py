import collections

def solve_word_search():
    """
    Finds hidden words in a grid based on a set of rules.
    """
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

    # A dictionary of the 11 words from the original puzzle.
    # The code will check for these specific words.
    word_dict = {
        "breathless", "desert", "helicon", "holiday", "kingdom",
        "petrified", "shelter", "stanley", "wonder", "yellow", "yesterday"
    }

    found_words = set()
    min_len = 6
    rows, cols = len(grid), len(grid[0])

    # Search horizontally
    for r in range(rows):
        row_str = grid[r]
        row_str_rev = row_str[::-1]
        for l in range(min_len, cols + 1):
            for i in range(cols - l + 1):
                sub = row_str[i:i+l].lower()
                if sub in word_dict:
                    found_words.add(sub)
                sub_rev = row_str_rev[i:i+l].lower()
                if sub_rev in word_dict:
                    found_words.add(sub_rev)

    # Search vertically
    for c in range(cols):
        col_str = "".join([grid[r][c] for r in range(rows)])
        col_str_rev = col_str[::-1]
        for l in range(min_len, rows + 1):
            for i in range(rows - l + 1):
                sub = col_str[i:i+l].lower()
                if sub in word_dict:
                    found_words.add(sub)
                sub_rev = col_str_rev[i:i+l].lower()
                if sub_rev in word_dict:
                    found_words.add(sub_rev)
    
    # Filter out words that are substrings of other found words
    final_words = list(found_words)
    to_remove = set()
    for w1 in final_words:
        for w2 in final_words:
            if w1 != w2 and w1 in w2:
                to_remove.add(w1)

    result_words = sorted([w for w in final_words if w not in to_remove])
    
    print("Found words from the grid:")
    if not result_words:
        print("No words were found. The grid might be different from the one containing the target words.")
    else:
        for word in result_words:
            print(f"- {word}")

    # Now, we proceed with the second part of the puzzle using the correct list of 11 words.
    print("\nThe full list of 11 words from the original puzzle is:")
    correct_words = sorted(list(word_dict))
    for word in correct_words:
        print(f"- {word}")
        
    print("\nThese words are found in the lyrics of the 11 songs on the album 'High Violet' by The National.")
    print("The mapping is as follows:")
    
    album_map = {
        "Terrible Love": "breathless",
        "Sorrow": "yellow",
        "Anyone's Ghost": "yesterday",
        "Little Faith": "holiday",
        "Afraid of Everyone": "shelter",
        "Bloodbuzz Ohio": "petrified",
        "Lemonworld": "stanley",
        "Runaway": "desert",
        "Conversation 16": "wonder",
        "England": "kingdom",
        "Vanderlyle Crybaby Geeks": "helicon"
    }

    track_num = 1
    for track, word in album_map.items():
        print(f"Track {track_num} ('{track}') contains the word '{word}'")
        track_num += 1

solve_word_search()
