from collections import deque

# A small set of words for demonstration purposes
# In a real scenario, this should be a comprehensive dictionary of English words
word_set = set(["VEND", "BEND", "BAND", "HAND", "HAST", "HEND", "HAST", "HANT", "HAST"])

def is_valid_word(word):
    return word in word_set

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            if c != word[i]:
                new_word = word[:i] + c + word[i+1:]
                if is_valid_word(new_word):
                    neighbors.append(new_word)
    return neighbors

def word_ladder(start, end):
    queue = deque([(start, [start])])
    visited = set([start])
    
    while queue:
        current_word, path = queue.popleft()
        
        if current_word == end:
            return path
        
        for neighbor in get_neighbors(current_word):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return []

# Find the word ladder
ladder = word_ladder("VEND", "HAST")
print(",".join(ladder))