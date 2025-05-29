# Final attempt with a different path
path = ['TOWS', 'TONS', 'TOIL', 'SOIL']

valid = True
for i in range(len(path)-1):
    if sum(1 for a, b in zip(path[i], path[i+1]) if a != b) != 1:
        valid = False
        break

print(','.join(path) if valid else "No valid path found")