from collections import deque

def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(position, occupied_positions):
    moves = [(2, 1), (2, -1), (-2, 1), (-2, -1),
             (1, 2), (1, -2), (-1, 2), (-1, -2)]
    possible_moves = []
    for dx, dy in moves:
        new_x = chr(ord(position[0]) + dx)
        new_y = str(int(position[1]) + dy)
        new_position = new_x + new_y
        if 'A' <= new_x <= 'D' and '1' <= new_y <= '3' and new_position not in occupied_positions:
            possible_moves.append(new_position)
    return possible_moves

def solve_knight_swap():
    initial_state = (('A1', 'B1'), ('B3', 'A2'))
    target_state = (('B3', 'A2'), ('A1', 'B1'))
    queue = deque([(initial_state, [], 'w')])
    visited = set()
    visited.add(initial_state)

    while queue:
        (white_knights, black_knights), moves, turn = queue.popleft()
        
        if (white_knights, black_knights) == target_state:
            return moves
        
        if turn == 'w':
            for i, knight in enumerate(white_knights):
                for move in get_possible_moves(knight, white_knights + black_knights):
                    new_white_knights = list(white_knights)
                    new_white_knights[i] = move
                    new_state = (tuple(new_white_knights), black_knights)
                    if new_state not in visited:
                        visited.add(new_state)
                        queue.append((new_state, moves + [f"w,{knight},{move}"], 'B'))
        else:
            for i, knight in enumerate(black_knights):
                for move in get_possible_moves(knight, white_knights + black_knights):
                    new_black_knights = list(black_knights)
                    new_black_knights[i] = move
                    new_state = (white_knights, tuple(new_black_knights))
                    if new_state not in visited:
                        visited.add(new_state)
                        queue.append((new_state, moves + [f"B,{knight},{move}"], 'w'))
    
    return "No"

# Execute the function and print the result
result = solve_knight_swap()
print(result)