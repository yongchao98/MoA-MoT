import itertools

def is_valid(p, n=5):
    k = len(p)
    rows = list(range(k))
    cols = p
    
    diffs = set()
    sums = set()
    
    for i in range(k):
        diffs.add((rows[i] - cols[i]) % n)
        sums.add((rows[i] + cols[i]) % n)
        
    return len(diffs) == k and len(sums) == k

count = 0
perms = list(itertools.permutations(range(4)))
valid_perms = []
for p in perms:
    if is_valid(p):
        count += 1
        valid_perms.append(p)

# print(f"Number of valid permutations for a fixed 4x4 subgrid: {count}")
# print(f"Valid permutations: {valid_perms}")