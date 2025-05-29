from collections import deque

def get_valid_words():
    # Common 4-letter English words
    words = set([
        'FOHN', 'GOBO', 'GOHN', 'GOON', 'GOBO', 'GOBS',
        'FOND', 'FONT', 'FOOT', 'FORD', 'FORK', 'FORM',
        'FORT', 'FOUL', 'FOUR', 'FOWL', 'FOXY', 'GOLD',
        'GONE', 'GONG', 'GOOD', 'GOOF', 'HOOK', 'HOOP',
        'HOPE', 'HORN', 'HOST', 'HOUR', 'HOWL', 'HOBO'
    ])
    return words

def differ_by_one(word1, word2):
    diff = 0
    for c1, c2 in zip(word1, word2):
        if c1 != c2:
            diff += 1
    return diff == 1

def find_path(start, end):
    words = get_valid_words()
    queue = deque([(start, [start])])
    seen = {start}
    
    while queue:
        word, path = queue.popleft()
        if word == end:
            return path
            
        for next_word in words:
            if next_word not in seen and differ_by_one(word, next_word):
                seen.add(next_word)
                queue.append((next_word, path + [next_word]))
    return None

path = find_path('FOHN', 'GOBO')
if path:
    print(','.join(path))
else:
    print("No valid path found")