# Let's try with a different set of intermediate words
word_list = {'NOSE', 'NOSS', 'SOSS', 'SUQS'}

def get_one_letter_changes(word):
    result = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            if c != word[i]:
                new_word = word[:i] + c + word[i+1:]
                if new_word in word_list:
                    result.append(new_word)
    return result

def find_path(start, target):
    if start not in word_list or target not in word_list:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current_word, path = queue.popleft()
        
        if current_word == target:
            return path
        
        for next_word in get_one_letter_changes(current_word):
            if next_word not in visited:
                visited.add(next_word)
                new_path = path + [next_word]
                queue.append((next_word, new_path))
    
    return None

from collections import deque
path = find_path('NOSE', 'SUQS')

if path:
    print(','.join(path))
else:
    print("No solution found")