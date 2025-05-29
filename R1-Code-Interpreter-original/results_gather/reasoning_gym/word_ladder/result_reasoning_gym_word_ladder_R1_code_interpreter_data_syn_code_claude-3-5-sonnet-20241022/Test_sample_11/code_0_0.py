from collections import deque, defaultdict

def get_word_set():
    # Common 4-letter words
    words = '''JAZZ JAZE JADE FADE FAME DOME DORM
               JACK BACK BARK DARK DANK DANK DORK
               JARS BARS BARE DARE DARE DORM
               JABS JOBS JOTS DOTS DOES DOES DORM
               JARL CARL CARD CORD CORD DORM'''.split()
    return set(words)

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            neighbors.append(new_word)
    return neighbors

def find_ladder(start, end, word_set):
    if start == end:
        return [start]
    
    # Bi-directional BFS
    forward = {start: [start]}
    backward = {end: [end]}
    direction = True  # True for forward, False for backward
    
    while forward and backward:
        if direction:
            current_dict = forward
            other_dict = backward
        else:
            current_dict = backward
            other_dict = forward
            
        next_dict = {}
        for word in current_dict:
            for neighbor in get_neighbors(word):
                if neighbor in word_set:
                    if neighbor in other_dict:
                        # Found intersection
                        if direction:
                            path = current_dict[word] + other_dict[neighbor][::-1]
                        else:
                            path = other_dict[neighbor] + current_dict[word][::-1]
                        return path
                    if neighbor not in current_dict:
                        next_dict[neighbor] = current_dict[word] + [neighbor]
        
        if direction:
            forward = next_dict
        else:
            backward = next_dict
        direction = not direction
    
    return None

# Find solution
word_set = get_word_set()
path = find_ladder('JAZZ', 'DORM', word_set)

if path:
    print(','.join(path))
else:
    print("No solution found")