from nltk.corpus import words

# Ensure the word list is available
word_list = set(words.words())

def is_valid_word(word):
    return word.lower() in word_list

def find_word_ladder(start, end):
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

ladder = find_word_ladder('RICH', 'NITE')
print(','.join(ladder))