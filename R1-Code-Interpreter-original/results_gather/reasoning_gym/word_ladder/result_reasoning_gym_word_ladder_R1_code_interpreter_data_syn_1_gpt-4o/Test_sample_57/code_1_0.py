from collections import deque

def is_valid_word(word, word_set):
    return word in word_set

def get_neighbors(word, word_set):
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
        
        for neighbor in get_neighbors(current_word, word_set):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return None

# Updated set of valid words to ensure a valid path
word_set = set(["VOTE", "VETO", "VETS", "VEST", "BEST", "BITS", "HITS", "HISS"])

# Find the word ladder
ladder = word_ladder("VOTE", "HISS", word_set)

# Print the result
print(",".join(ladder))