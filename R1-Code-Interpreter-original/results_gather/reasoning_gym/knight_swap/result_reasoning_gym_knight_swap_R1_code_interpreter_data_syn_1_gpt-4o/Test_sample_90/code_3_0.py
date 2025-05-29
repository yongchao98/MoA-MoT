from collections import deque

def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def swap_knights():
    # Initial positions
    white_knights = ['C1', 'D1']
    black_knights = ['A3', 'A2']
    empty_squares = ['A1', 'B1', 'B2', 'B3', 'B4', 'C2', 'C3', 'C4', 'D2', 'D3']
    
    # Target positions
    target_white = ['A3', 'A2']
    target_black = ['C1', 'D1']
    
    # Possible moves for a knight
    knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]
    
    # Convert position to coordinates
    def pos_to_coord(pos):
        return (ord(pos[0]) - ord('A'), int(pos[1]) - 1)
    
    # Convert coordinates to position
    def coord_to_pos(coord):
        return chr(coord[0] + ord('A')) + str(coord[1] + 1)
    
    # BFS queue
    queue = deque()
    queue.append((white_knights, black_knights, empty_squares, []))
    
    # BFS to find the solution
    while queue:
        w_knights, b_knights, empties, moves = queue.popleft()
        
        # Check if the swap is complete
        if set(w_knights) == set(target_white) and set(b_knights) == set(target_black):
            return moves
        
        # Generate all possible moves for white knights
        for i, w_knight in enumerate(w_knights):
            w_coord = pos_to_coord(w_knight)
            for move in knight_moves:
                new_w_coord = (w_coord[0] + move[0], w_coord[1] + move[1])
                new_w_pos = coord_to_pos(new_w_coord)
                if new_w_pos in empties:
                    new_w_knights = w_knights[:i] + [new_w_pos] + w_knights[i+1:]
                    new_empties = [w_knight if e == new_w_pos else e for e in empties]
                    queue.append((new_w_knights, b_knights, new_empties, moves + [f"w,{w_knight},{new_w_pos}"]))
        
        # Generate all possible moves for black knights
        for i, b_knight in enumerate(b_knights):
            b_coord = pos_to_coord(b_knight)
            for move in knight_moves:
                new_b_coord = (b_coord[0] + move[0], b_coord[1] + move[1])
                new_b_pos = coord_to_pos(new_b_coord)
                if new_b_pos in empties:
                    new_b_knights = b_knights[:i] + [new_b_pos] + b_knights[i+1:]
                    new_empties = [b_knight if e == new_b_pos else e for e in empties]
                    queue.append((w_knights, new_b_knights, new_empties, moves + [f"B,{b_knight},{new_b_pos}"]))
    
    return "No"

# Execute the function and print the result
print(swap_knights())