def is_every_crop_hydrated(field):
    rows = len(field)
    cols = len(field[0])
    
    # Create a hydrated map initialized to False
    hydrated = [[False] * cols for _ in range(rows)]
    
    # Directions for the 8 surrounding tiles
    directions = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    
    # Mark hydrated areas around each water source
    for r in range(rows):
        for c in range(cols):
            if field[r][c] == 'w':
                # Mark the water source itself as hydrated
                hydrated[r][c] = True
                # Mark the surrounding 8 tiles
                for dr, dc in directions:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        hydrated[nr][nc] = True
    
    # Check if all crops are hydrated
    for r in range(rows):
        for c in range(cols):
            if field[r][c] == 'c' and not hydrated[r][c]:
                return False
    
    return True

# Given input
field = [['c', 'c', 'c', 'c', 'c', 'c'], 
         ['c', 'c', 'w', 'c', 'c', 'c'], 
         ['c', 'c', 'c', 'w', 'c', 'c'], 
         ['c', 'c', 'c', 'c', 'c', 'c'], 
         ['c', 'w', 'c', 'c', 'c', 'w']]

# Determine if every crop is hydrated
result = is_every_crop_hydrated(field)
print(result)