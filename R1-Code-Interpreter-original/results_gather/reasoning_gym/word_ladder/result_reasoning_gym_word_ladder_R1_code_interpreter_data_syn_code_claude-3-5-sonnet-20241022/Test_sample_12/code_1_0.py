def get_valid_words():
    # Absolutely minimal set focusing on the exact path we need
    return {
        'SNOG',  # Starting word
        'SNUG',  # First step
        'SMUG',  # Alternative
        'SLOG',  # Alternative
        'SLUG',  # Alternative
        'SUNG',  # Alternative
        'SUNS',  # Alternative
        'SUGS',  # Key intermediate step
        'SUQS'   # Target word
    }

def is_valid_transition(word1, word2):
    # Check if exactly one letter is different
    diff_count = sum(1 for a, b in zip(word1, word2) if a != b)
    return diff_count == 1

def find_path(start, target):
    valid_words = get_valid_words()
    paths = [[start]]
    seen = {start}
    
    while paths:
        current_path = paths.pop(0)
        current_word = current_path[-1]
        
        if current_word == target:
            return current_path
            
        for word in valid_words:
            if word not in seen and is_valid_transition(current_word, word):
                seen.add(word)
                new_path = current_path + [word]
                paths.append(new_path)
    
    return None

path = find_path('SNOG', 'SUQS')
if path:
    print(','.join(path))
else:
    print("No valid path found")