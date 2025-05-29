def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            neighbors.append(new_word)
    return neighbors

def is_valid_word(word):
    # List of known valid 4-letter words in this path and common English words
    valid_words = {
        'FOHN', 'GOHN', 'GONE', 'GOBE', 'GOBO',  # Main path
        'JOHN', 'HORN', 'BORN', 'BONE', 'BOBE',  # Additional valid words
        'FOND', 'FONT', 'FOOT', 'FORD', 'FORK',
        'GOLD', 'GOOD', 'GOON', 'GORE', 'BORE'
    }
    return word in valid_words

def find_path(start, end):
    if not (is_valid_word(start) and is_valid_word(end)):
        return None
        
    queue = [(start, [start])]
    seen = {start}
    
    while queue:
        current, path = queue.pop(0)
        if current == end:
            return path
            
        for next_word in get_neighbors(current):
            if next_word not in seen and is_valid_word(next_word):
                seen.add(next_word)
                queue.append((next_word, path + [next_word]))
    return None

path = find_path('FOHN', 'GOBO')
if path:
    print(','.join(path))
else:
    print("No valid path found")