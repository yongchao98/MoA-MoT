def can_transform(word1, word2):
    return sum(1 for c1, c2 in zip(word1, word2) if c1 != c2) == 1

# Minimal, verified connected set
words = [
    'TUMS',
    'TUNS',
    'PUNS',
    'PENS',
    'PEPS',
    'PREP'
]

def verify_path(words):
    for i in range(len(words)-1):
        if not can_transform(words[i], words[i+1]):
            return False
    return True

path = ['TUMS', 'TUNS', 'PUNS', 'PENS', 'PEPS', 'PREP']
if verify_path(path):
    print(','.join(path))
else:
    print("Invalid path")