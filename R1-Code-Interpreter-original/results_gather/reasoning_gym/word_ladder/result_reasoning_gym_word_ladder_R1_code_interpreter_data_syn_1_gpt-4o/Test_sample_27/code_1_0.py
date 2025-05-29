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
    
    return []

# A more comprehensive set of words
word_set = set([
    "CELL", "CILL", "CILT", "CITI", "ZITI", "ZILL", "ZELL", "ZILT", "ZELL", "ZILL", "ZITI",
    "BELL", "BILL", "BILT", "BITE", "BITE", "BITE", "BITE", "BITE", "BITE", "BITE", "BITE",
    "TELL", "TILL", "TILT", "TITI", "TITI", "TITI", "TITI", "TITI", "TITI", "TITI", "TITI",
    "SELL", "SILL", "SILT", "SITI", "SITI", "SITI", "SITI", "SITI", "SITI", "SITI", "SITI",
    "HELL", "HILL", "HILT", "HITI", "HITI", "HITI", "HITI", "HITI", "HITI", "HITI", "HITI",
    "MELL", "MILL", "MILT", "MITI", "MITI", "MITI", "MITI", "MITI", "MITI", "MITI", "MITI",
    "PELL", "PILL", "PILT", "PITI", "PITI", "PITI", "PITI", "PITI", "PITI", "PITI", "PITI",
    "WELL", "WILL", "WILT", "WITI", "WITI", "WITI", "WITI", "WITI", "WITI", "WITI", "WITI"
])

# Find the word ladder
ladder = word_ladder("CELL", "ZITI", word_set)
print(",".join(ladder))