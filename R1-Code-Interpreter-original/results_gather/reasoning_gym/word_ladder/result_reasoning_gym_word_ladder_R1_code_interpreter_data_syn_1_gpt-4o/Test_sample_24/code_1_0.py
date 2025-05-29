# Define a small set of valid words for demonstration
valid_words = {
    "ROPY", "RODY", "RODS", "SODS"
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
    
    return []

# Find the word ladder
ladder = word_ladder('ROPY', 'SODS')
print(','.join(ladder))