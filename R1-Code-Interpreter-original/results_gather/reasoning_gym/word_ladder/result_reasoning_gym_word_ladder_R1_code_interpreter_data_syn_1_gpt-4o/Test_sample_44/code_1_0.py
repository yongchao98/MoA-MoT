# Define a small set of valid words for the transformation
valid_words = {
    "RICH", "RICE", "RITE", "RIME", "NICE", "NITE", "NICE", "NITE", "RICH", "RICE", "RIME", "RITE", "NICE", "NITE"
}

def is_valid_word(word):
    return word in valid_words

def word_ladder(start, end):
    from collections import deque
    
    # Queue for BFS
    queue = deque([(start, [start])])
    
    while queue:
        current_word, path = queue.popleft()
        
        if current_word == end:
            return path
        
        for i in range(len(current_word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if c != current_word[i]:
                    next_word = current_word[:i] + c + current_word[i+1:]
                    if is_valid_word(next_word) and next_word not in path:
                        queue.append((next_word, path + [next_word]))
    
    return None

# Find the word ladder
ladder = word_ladder('RICH', 'NITE')
print(','.join(ladder))