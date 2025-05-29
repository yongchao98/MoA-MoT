from collections import deque

def is_valid_word(word, word_list):
    return word in word_list

def word_ladder(start, end, word_list):
    queue = deque([(start, [start])])
    word_set = set(word_list)
    
    while queue:
        current_word, path = queue.popleft()
        
        if current_word == end:
            return path
        
        for i in range(len(current_word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                next_word = current_word[:i] + c + current_word[i+1:]
                if next_word in word_set and next_word != current_word:
                    queue.append((next_word, path + [next_word]))
                    word_set.remove(next_word)
    
    return []

# List of valid English words (for demonstration purposes, a small subset)
word_list = ["USED", "USER", "CURE", "CUES", "CUES", "CUBE", "CUBE", "CUTE", "CITE", "SITE", "SILE", "FILE"]

# Find the word ladder
ladder = word_ladder("USED", "FILE", word_list)
print(",".join(ladder))