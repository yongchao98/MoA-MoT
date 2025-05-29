from collections import deque

# A small set of valid words for demonstration purposes
valid_words = set(["LOGS", "GOGS", "GOES", "GORE", "LORE", "CORE", "BORE", "MORE", "SORE", "FORE"])

def is_valid_word(word):
    return word in valid_words

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            if c != word[i]:
                new_word = word[:i] + c + word[i+1:]
                if is_valid_word(new_word):
                    neighbors.append(new_word)
    return neighbors

def find_word_ladder(start, end):
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
    
    return None

ladder = find_word_ladder("LOGS", "GORE")
if ladder:
    print(",".join(ladder))
else:
    print("No valid word ladder found.")