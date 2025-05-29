from collections import deque

def get_valid_words():
    # Common 4-letter words
    words = set()
    try:
        with open('/usr/share/dict/words', 'r') as f:
            for line in f:
                word = line.strip().upper()
                if len(word) == 4:
                    words.add(word)
    except:
        # Fallback to a smaller set of common 4-letter words
        words = {'SNOG', 'SNUG', 'SLUG', 'SLUM', 'SLUR', 'SPUR', 'SPUN', 'SPIN',
                'SUNG', 'SUNS', 'SUQS', 'SUDS', 'SUBS', 'SUES', 'SUNG'}
    return words

def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            neighbors.append(new_word)
    return neighbors

def find_path(start, end, valid_words):
    queue = deque([(start, [start])])
    seen = {start}
    
    while queue:
        word, path = queue.popleft()
        if word == end:
            return path
            
        for neighbor in get_neighbors(word):
            if neighbor in valid_words and neighbor not in seen:
                seen.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    return None

valid_words = get_valid_words()
path = find_path('SNOG', 'SUQS', valid_words)
if path:
    print(','.join(path))
else:
    print("No valid path found")