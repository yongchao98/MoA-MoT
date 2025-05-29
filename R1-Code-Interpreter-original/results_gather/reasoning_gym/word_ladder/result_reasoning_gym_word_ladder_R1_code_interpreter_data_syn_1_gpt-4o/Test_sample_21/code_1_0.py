from collections import deque

def is_valid_word(word, word_set):
    return word in word_set

def generate_neighbors(word, word_set):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            if c != word[i]:
                new_word = word[:i] + c + word[i+1:]
                if is_valid_word(new_word, word_set):
                    neighbors.append(new_word)
    return neighbors

def word_ladder(start, end, word_set):
    queue = deque([(start, [start])])
    visited = set([start])
    
    while queue:
        current_word, path = queue.popleft()
        
        if current_word == end:
            return path
        
        for neighbor in generate_neighbors(current_word, word_set):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return []

# A more comprehensive set of words for demonstration purposes
word_set = set([
    "READ", "ROAD", "ROOD", "ROOF", "ROOK", "ROCK", "SOCK", "SICK", "SHIP", 
    "SHED", "SHAD", "SHID", "SHAP", "SLIP", "SNIP", "SNAP", "SNAG", "SHAG", 
    "SHIP", "SEED", "SEEP", "SHEP", "SHIP"
])

# Find the word ladder
ladder = word_ladder("READ", "SHIP", word_set)
print(",".join(ladder))