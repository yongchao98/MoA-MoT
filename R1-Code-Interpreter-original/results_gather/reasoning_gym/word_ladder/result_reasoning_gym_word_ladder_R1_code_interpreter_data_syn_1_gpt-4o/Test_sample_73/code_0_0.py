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

def word_ladder(start, end, word_list):
    word_set = set(word_list)
    queue = deque([(start, [start])])
    while queue:
        current_word, path = queue.popleft()
        if current_word == end:
            return path
        for neighbor in get_neighbors(current_word, word_set):
            if neighbor not in path:
                queue.append((neighbor, path + [neighbor]))
    return []

# Example word list (a small subset of valid words for demonstration)
word_list = ["TREE", "FREE", "FLEE", "FLEW", "FLOW", "CLOW", "CLOWN", "COWS"]

# Find the word ladder
ladder = word_ladder("TREE", "COWS", word_list)
print(",".join(ladder))