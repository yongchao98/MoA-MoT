from collections import deque

# A more comprehensive set of words for demonstration purposes
word_set = set([
    "DIBS", "DIMS", "DUMP", "LUMP", "LIMB", "LIMPS", "LIMN", "LUMP", 
    "DUMB", "DAMP", "LAMP", "LIMP", "LUMP"
])

def is_valid_word(word):
    return word in word_set

def word_ladder(start, end):
    queue = deque([(start, [start])])
    while queue:
        current_word, path = queue.popleft()
        if current_word == end:
            return path
        for i in range(len(current_word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                next_word = current_word[:i] + c + current_word[i+1:]
                if is_valid_word(next_word) and next_word not in path:
                    queue.append((next_word, path + [next_word]))
    return []

# Find the word ladder
ladder = word_ladder("DIBS", "LUMP")
print(",".join(ladder))