def find_swapped_lines():
    """
    Analyzes the poem to find the two lines that have been swapped,
    breaking the sestina form.
    """
    
    # The six key end-words from the first stanza (lines 1-6)
    # 1. mountains, 2. valleys, 3. forests, 4. music, 5. morning, 6. evening
    end_words_stanza_1_to_6 = [
        # Stanza 1 (Correct by definition)
        ['mountains', 'valleys', 'forests', 'music', 'morning', 'evening'],
        # Stanza 2 (Expected: 6-1-5-2-4-3)
        ['evening', 'mountains', 'morning', 'valleys', 'music', 'forests'],
        # Stanza 3 (Expected: 3-6-4-1-2-5)
        ['forests', 'evening', 'music', 'mountains', 'valleys', 'morning'],
        # Stanza 4 (Expected: 5-3-2-6-1-4)
        ['morning', 'forests', 'valleys', 'evening', 'mountains', 'music'],
        # Stanza 5 (Expected: 4-5-1-3-6-2)
        ['music', 'morning', 'mountains', 'forests', 'evening', 'valleys'],
        # Stanza 6 (Expected: 2-4-6-5-3-1)
        ['valleys', 'music', 'evening', 'morning', 'forests', 'mountains']
    ]

    # The expected pattern repeats for the second half (a double sestina)
    end_words_stanza_7_to_12 = end_words_stanza_1_to_6

    all_expected_patterns = end_words_stanza_1_to_6 + end_words_stanza_7_to_12

    # The actual end-words from the provided poem text
    actual_end_words = [
        "mountains", "valleys", "forests", "music", "morning", "evening",      # Stanza 1 (1-6)
        "evening", "mountains", "morning", "valleys", "music", "forests",      # Stanza 2 (7-12)
        "forests", "evening", "music", "mountains", "valleys", "morning",      # Stanza 3 (13-18)
        "morning", "forests", "valleys", "evening", "mountains", "music",      # Stanza 4 (19-24)
        "music", "forests", "mountains", "morning", "evening", "valleys",      # Stanza 5 (25-30) -> Contains the error
        "valleys", "music", "evening", "morning", "forests", "mountains",      # Stanza 6 (31-36)
        "mountains", "valleys", "forests", "music", "morning", "evening",      # Stanza 7 (37-42)
        "evening", "mountains", "morning", "valleys", "music", "forests",      # Stanza 8 (43-48)
        "forests", "evening", "music", "mountains", "valleys", "morning",      # Stanza 9 (49-54)
        "morning", "forests", "valleys", "evening", "mountains", "music",      # Stanza 10 (55-60)
        "music", "morning", "mountains", "forests", "evening", "valleys",      # Stanza 11 (61-66)
        "valleys", "music", "evening", "morning", "forests", "mountains"       # Stanza 12 (67-72)
    ]
    
    # Find the stanza with the discrepancy
    for i in range(12):
        expected = all_expected_patterns[i]
        actual = actual_end_words[i*6 : (i+1)*6]
        
        if expected != actual:
            stanza_num = i + 1
            start_line = i * 6 + 1
            
            # Find the indices of the mismatched words within the stanza
            mismatched_indices = []
            for j in range(6):
                if expected[j] != actual[j]:
                    mismatched_indices.append(j)
            
            # The line numbers are the start line + the 0-based index
            line1 = start_line + mismatched_indices[0]
            line2 = start_line + mismatched_indices[1]
            
            print(f"{line1} and {line2}")
            return

find_swapped_lines()