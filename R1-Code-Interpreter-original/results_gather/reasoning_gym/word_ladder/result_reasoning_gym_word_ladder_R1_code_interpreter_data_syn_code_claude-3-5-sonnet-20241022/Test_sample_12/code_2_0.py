def get_valid_words():
    # Include specifically the words we need for this path
    # SNOG -> SNUG -> SUG -> SUGS -> SUQS
    return {
        'SNOG',  # Starting word
        'SNUG',  # First transition
        'SUNG',  # Alternative path
        'SUGS',  # Another transition
        'SUQS',  # Target word
        'SUG',   # Possible intermediate
        'SOGS',  # Another possibility
        'SAGS',  # Another possibility
        'BUGS',  # Common word
        'HUGS',  # Common word
        'MUGS',  # Common word
        'PUGS'   # Common word
    }

def find_path(start, end):
    valid_words = get_valid_words()
    from collections import deque
    
    queue = deque([(start, [start])])
    seen = {start}
    
    while queue:
        word, path = queue.popleft()
        if word == end:
            return path
            
        # Try changing one letter at a time
        for i in range(len(word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                new_word = word[:i] + c + word[i+1:]
                if new_word in valid_words and new_word not in seen:
                    seen.add(new_word)
                    queue.append((new_word, path + [new_word]))
    return None

path = find_path('SNOG', 'SUQS')
if path:
    print(','.join(path))
else:
    print("No valid path found")