# Test final verified path
path = ['SHIV', 'SHIP', 'SLIP', 'SLAP', 'SNAP', 'TRAP', 'TARP']

def verify_path(path):
    for i in range(len(path)-1):
        diff = sum(1 for a, b in zip(path[i], path[i+1]) if a != b)
        if diff != 1:
            return False
    return True

if verify_path(path):
    print(','.join(path))
else:
    print("No valid path found")