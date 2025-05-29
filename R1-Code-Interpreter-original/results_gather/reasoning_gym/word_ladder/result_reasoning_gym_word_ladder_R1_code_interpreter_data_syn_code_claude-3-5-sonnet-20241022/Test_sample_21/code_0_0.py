from collections import deque, defaultdict

# Define a function to check if two words differ by exactly one letter
def diff_by_one(word1, word2):
    return sum(1 for i in range(len(word1)) if word1[i] != word2[i]) == 1

# Create a set of common 4-letter words
words = {
    'READ', 'ROAD', 'RAID', 'RAIL', 'SAIL', 'SHIP',
    'HEAD', 'HEAL', 'HEAR', 'HEAT', 'HEAP',
    'LEAD', 'LOAD', 'LEAN', 'LEAP',
    'REAL', 'SEAL', 'SEAT', 'SLIP', 'SHOP',
    'REED', 'SEED', 'SEND', 'SHED', 'SHAD',
    'REAP', 'HEAP', 'HELP', 'HELD',
    'READ', 'ROAD', 'ROAM', 'FOAM', 'FEAR',
}

# Build graph of word connections
graph = defaultdict(list)
for word1 in words:
    for word2 in words:
        if diff_by_one(word1, word2):
            graph[word1].append(word2)

# BFS to find shortest path
def find_path(start, end):
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        word = path[-1]
        
        if word == end:
            return path
            
        for next_word in graph[word]:
            if next_word not in visited:
                visited.add(next_word)
                queue.append(path + [next_word])
    
    return None

# Find and print the path
path = find_path('READ', 'SHIP')
if path:
    print(','.join(path))
else:
    print("No path found")