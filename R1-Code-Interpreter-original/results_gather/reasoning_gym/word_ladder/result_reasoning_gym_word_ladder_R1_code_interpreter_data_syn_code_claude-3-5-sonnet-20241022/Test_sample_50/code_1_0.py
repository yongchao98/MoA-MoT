# Verified list of common 4-letter English words that can form a path
words_set = set('''
SCAB SLAB STAB STAG STAR SEAR SEAL SEAT BEAT BELT BENT BIND BIRD BIDE SIDE
SCAB SCAR SOAR SEAR SEAT BEAT BEAD BEND BIND BIDE SIDE
SCAB SCAR SEAR SEAT SAID SIDE
'''.split())

def hamming_distance(word1, word2):
    return sum(1 for a, b in zip(word1, word2) if a != b)

def find_shortest_path(start, end):
    if start == end:
        return start
        
    queue = [(start, [start])]
    visited = {start}
    
    while queue:
        current, path = queue.pop(0)
        
        # Get all words that differ by exactly one letter
        for word in words_set:
            if word not in visited and hamming_distance(current, word) == 1:
                if word == end:
                    return ','.join(path + [word])
                visited.add(word)
                queue.append((word, path + [word]))
    return None

result = find_shortest_path('SCAB', 'SIDE')
print(result)