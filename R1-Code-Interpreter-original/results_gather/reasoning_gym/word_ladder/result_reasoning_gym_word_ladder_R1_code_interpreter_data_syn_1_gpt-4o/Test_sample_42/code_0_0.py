from collections import deque

def is_valid_word(word, word_set):
    return word in word_set

def word_ladder(start, end, word_list):
    word_set = set(word_list)
    queue = deque([(start, [start])])
    
    while queue:
        current_word, path = queue.popleft()
        
        if current_word == end:
            return path
        
        for i in range(len(current_word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                next_word = current_word[:i] + c + current_word[i+1:]
                if next_word != current_word and is_valid_word(next_word, word_set):
                    queue.append((next_word, path + [next_word]))
                    word_set.remove(next_word)  # Remove to prevent revisiting
    
    return []

# Define a list of valid English words (for demonstration purposes, a small subset)
valid_words = ["LESS", "LENS", "LEND", "LAND", "DAND", "DANG"]

# Find the word ladder
ladder = word_ladder("LESS", "DANG", valid_words)
print(",".join(ladder))