def is_valid_word(word, word_list):
    return word in word_list

def word_ladder(start, end, word_list):
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
                    if is_valid_word(next_word, word_list) and next_word not in path:
                        queue.append((next_word, path + [next_word]))
    
    return None

# Define a small set of valid words for this problem
word_list = {"ICED", "ACED", "AGED", "AXED", "ARCH", "ACHE", "ARCS", "RICH", "RICE", "RIDE", "RIDE", "RICH"}

# Find the word ladder
ladder = word_ladder('ICED', 'ARCH', word_list)
print(ladder)