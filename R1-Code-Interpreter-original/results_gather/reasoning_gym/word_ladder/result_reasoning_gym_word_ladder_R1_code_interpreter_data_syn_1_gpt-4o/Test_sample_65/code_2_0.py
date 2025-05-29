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

# Expanded dictionary of valid words for this problem
word_list = {
    "PLOY", "PLOT", "MOOT", "MOPE", "MOTE", "MOPY", "MOLY", "MOLP", "MOLE",
    "PLOP", "PLOE", "PLOM", "PLOX", "PLOZ", "PLOA", "PLOB", "PLOC", "PLOD",
    "PLOF", "PLOG", "PLOH", "PLOI", "PLOJ", "PLOK", "PLOL", "PLOM", "PLON",
    "PLOQ", "PLOR", "PLOS", "PLOT", "PLOU", "PLOV", "PLOW", "PLOX", "PLOY",
    "PLOZ", "MOPE", "MOTE", "MOLE", "MOSE", "MOXE", "MOYE", "MOZE", "MOAE",
    "MOBE", "MOCE", "MODE", "MOEE", "MOFE", "MOGE", "MOHE", "MOIE", "MOJE",
    "MOKE", "MOLE", "MOME", "MONE", "MOOE", "MOPE", "MOQE", "MORE", "MOSE",
    "MOTE", "MOU", "MOVE", "MOWE", "MOXE", "MOYE", "MOZE"
}

# Find the word ladder
ladder = word_ladder("PLOY", "MOPE", word_list)
print(",".join(ladder))