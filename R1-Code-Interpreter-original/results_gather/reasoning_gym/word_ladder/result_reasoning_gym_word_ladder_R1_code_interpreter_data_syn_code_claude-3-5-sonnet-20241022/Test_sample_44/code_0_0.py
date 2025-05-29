from collections import deque, defaultdict
import re

def get_word_list():
    # Get a list of 4-letter words from a basic English word list
    word_list = set()
    with open('/usr/share/dict/words', 'r') as f:
        for line in f:
            word = line.strip().upper()
            if len(word) == 4 and re.match('^[A-Z]+$', word):
                word_list.add(word)
    return word_list

def get_neighbors(word):
    # Generate all possible one-letter transformations
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            neighbors.append(new_word)
    return neighbors

def find_word_ladder(start, end, word_list):
    if start not in word_list or end not in word_list:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current_word, path = queue.popleft()
        if current_word == end:
            return path
        
        for neighbor in get_neighbors(current_word):
            if neighbor in word_list and neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return None

# Main execution
start_word = 'RICH'
end_word = 'NITE'
word_list = get_word_list()
path = find_word_ladder(start_word, end_word, word_list)

if path:
    print(','.join(path))
else:
    print("No valid path found")