def is_valid_word(word, valid_words):
    return word in valid_words

def transform_word_ladder(start, end, valid_words):
    from collections import deque
    
    # Initialize the queue with the start word
    queue = deque([(start, [start])])
    
    while queue:
        current_word, path = queue.popleft()
        
        if current_word == end:
            return path
        
        for i in range(len(current_word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if c != current_word[i]:
                    next_word = current_word[:i] + c + current_word[i+1:]
                    if is_valid_word(next_word, valid_words) and next_word not in path:
                        queue.append((next_word, path + [next_word]))
    
    return []

# Define a small set of valid words
valid_words = {
    'PROD', 'PROS', 'PODS', 'PARD', 'PARK', 'MARK', 'MARD', 'MARL'
}

# Find the transformation path
path = transform_word_ladder('PROD', 'MARL', valid_words)
print(','.join(path))