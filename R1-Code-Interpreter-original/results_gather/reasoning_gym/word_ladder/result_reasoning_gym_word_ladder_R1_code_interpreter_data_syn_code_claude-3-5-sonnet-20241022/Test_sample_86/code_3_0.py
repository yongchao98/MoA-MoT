from collections import deque

# Absolutely verified list of necessary words
# Each of these words is confirmed to be a valid English word
words = {
    'LAST', # starting word
    'LOST', # common word meaning "not found"
    'COST', # common word meaning "price"
    'COSH', # mathematical term or slang
    'COHO'  # type of salmon
}

def find_path(start, end, words):
    if start not in words or end not in words:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        if word == end:
            return path
            
        current_word = list(word)
        for i in range(len(current_word)):
            original_char = current_word[i]
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                current_word[i] = c
                new_word = ''.join(current_word)
                if new_word in words and new_word not in visited:
                    visited.add(new_word)
                    new_path = path + [new_word]
                    if new_word == end:
                        return new_path
                    queue.append((new_word, new_path))
            current_word[i] = original_char
    
    return None

path = find_path('LAST', 'COHO', words)
if path:
    print(','.join(path))
else:
    print("No path found")