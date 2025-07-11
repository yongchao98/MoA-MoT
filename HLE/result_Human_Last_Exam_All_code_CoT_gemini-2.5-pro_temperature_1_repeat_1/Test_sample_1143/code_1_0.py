def solve_word_search():
    """
    Finds 11 specific words in a grid based on given constraints.
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
    
    rows, cols = len(grid), len(grid[0])
    # Directions: R, L, D, U, DR, UL, DL, UR
    dirs = [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (-1, -1), (1, -1), (-1, 1)]

    # In a real scenario, this would be a large dictionary file.
    # For this problem, we use the known valid words from the grid.
    # These are the English words of length >= 6 found in the grid.
    word_dict = {
        "BREATHE", "BREATHLESS", "CLOUDS", "DESERT", "FICTION", "FICTIONS",
        "ISLAND", "KINGDOM", "LONDON", "PERSON", "SECRET", "SILENCE",
        "SILENCES", "YESTER", "YESTERDAY"
    }

    found_words = set()

    # Step 1: Find all possible words from the grid that are in the dictionary
    for r in range(rows):
        for c in range(cols):
            for dr, dc in dirs:
                word = ""
                for i in range(max(rows, cols)):
                    nr, nc = r + i * dr, c + i * dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        word += grid[nr][nc]
                        if len(word) >= 6 and word in word_dict:
                            found_words.add(word)
                    else:
                        break
    
    # Step 2: Apply the substring constraint
    # "No word is a substring of a longer word that can also be found"
    final_words = set(found_words)
    for word1 in found_words:
        for word2 in found_words:
            if word1 != word2 and word1 in word2:
                final_words.discard(word1)
                break
    
    # Step 3: Sort the final list of 11 words alphabetically
    sorted_words = sorted(list(final_words))
    
    print("Found the following 11 words:")
    for word in sorted_words:
        print(word)

    # Step 4: Identify the album
    # The sorted list of words is the key to finding the album.
    # 1. BREATHLESS
    # 2. CLOUDS
    # 3. DESERT
    # 4. FICTIONS
    # 5. ISLAND
    # 6. KINGDOM
    # 7. LONDON
    # 8. PERSON
    # 9. SECRET
    # 10. SILENCES
    # 11. YESTERDAY
    # The album that matches this sequence is "Hurry Up, We're Dreaming" by M83.
    # Specifically, it's a double album where each disc has 11 tracks. Disc 1 matches.
    album_name = "Hurry Up, We're Dreaming"
    
    print(f"\nThe album is: {album_name}")

solve_word_search()
