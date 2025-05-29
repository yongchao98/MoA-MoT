def is_valid_solution(positions):
    # Check if any queens can attack each other
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            r1, c1 = positions[i]
            r2, c2 = positions[j]
            
            # Check same row or column
            if r1 == r2 or c1 == c2:
                return False
                
            # Check diagonals
            if abs(r1 - r2) == abs(c1 - c2):
                return False
    
    return True

# Combine fixed positions with code output
positions = [
    (0, 1),  # fixed
    (3, 6),  # fixed
    (5, 7),  # fixed
    (0, 0),  # from code
    (1, 2),  # from code
    (2, 4),  # from code
    (3, 1),  # from code
    (4, 3)   # from code
]

# Sort positions for consistent output
positions.sort()

# Check if solution is valid and print result
if is_valid_solution(positions):
    result = ", ".join(f"{pos[0]} {pos[1]}" for pos in positions)
    print(f"<<{result}>>")
else:
    print("Invalid solution - queens can attack each other")