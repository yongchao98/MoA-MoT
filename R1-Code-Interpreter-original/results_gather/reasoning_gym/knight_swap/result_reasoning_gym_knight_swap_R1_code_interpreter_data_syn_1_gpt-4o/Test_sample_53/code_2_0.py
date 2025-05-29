from collections import deque

def is_valid_move(start, end, occupied):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx, dy) in [(2, 1), (1, 2)] and end not in occupied

def swap_knights():
    # Initial positions
    initial_white = ('B1', 'A2')
    initial_black = ('C3', 'B2')
    initial_empty = ('B3', 'C1')
    
    # Target positions
    target_white = ('C3', 'B2')
    target_black = ('B1', 'A2')
    
    # Queue for BFS
    queue = deque()
    queue.append((initial_white, initial_black, initial_empty, []))
    
    # Visited set to avoid cycles
    visited = set()
    visited.add((initial_white, initial_black))
    
    while queue:
        white_knights, black_knights, empty_squares, moves = queue.popleft()
        
        # Check if the target configuration is reached
        if white_knights == target_white and black_knights == target_black:
            return moves
        
        # Generate possible moves for white knights
        for i, white in enumerate(white_knights):
            for empty in empty_squares:
                if is_valid_move(white, empty, white_knights + black_knights):
                    new_white_knights = list(white_knights)
                    new_white_knights[i] = empty
                    new_empty_squares = list(empty_squares)
                    new_empty_squares.remove(empty)
                    new_empty_squares.append(white)
                    new_state = (tuple(new_white_knights), black_knights, tuple(new_empty_squares))
                    if new_state[:2] not in visited:
                        visited.add(new_state[:2])
                        queue.append((tuple(new_white_knights), black_knights, tuple(new_empty_squares), moves + [f"w,{white},{empty}"]))
        
        # Generate possible moves for black knights
        for i, black in enumerate(black_knights):
            for empty in empty_squares:
                if is_valid_move(black, empty, white_knights + black_knights):
                    new_black_knights = list(black_knights)
                    new_black_knights[i] = empty
                    new_empty_squares = list(empty_squares)
                    new_empty_squares.remove(empty)
                    new_empty_squares.append(black)
                    new_state = (white_knights, tuple(new_black_knights), tuple(new_empty_squares))
                    if new_state[:2] not in visited:
                        visited.add(new_state[:2])
                        queue.append((white_knights, tuple(new_black_knights), tuple(new_empty_squares), moves + [f"B,{black},{empty}"]))
    
    return "No"

# Execute the function and print the result
print(swap_knights())