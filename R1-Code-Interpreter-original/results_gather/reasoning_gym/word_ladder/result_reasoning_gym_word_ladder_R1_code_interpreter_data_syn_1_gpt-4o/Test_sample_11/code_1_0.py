from collections import deque

def is_valid_word(word, word_list):
    return word in word_list

def get_neighbors(word, word_list):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            if c != word[i]:
                new_word = word[:i] + c + word[i+1:]
                if is_valid_word(new_word, word_list):
                    neighbors.append(new_word)
    return neighbors

def word_ladder(start, end, word_list):
    queue = deque([(start, [start])])
    visited = set([start])
    
    while queue:
        current_word, path = queue.popleft()
        
        if current_word == end:
            return path
        
        for neighbor in get_neighbors(current_word, word_list):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return []

# A small dictionary of valid words for demonstration purposes
word_list = set(["JAZZ", "JARS", "JAMS", "JUMP", "DUMP", "DAMP", "DORM", "FORM", "FARM", "FIRM", "WORM", "WARM", "WORD"])

# Find the word ladder
ladder = word_ladder("JAZZ", "DORM", word_list)
print(",".join(ladder))