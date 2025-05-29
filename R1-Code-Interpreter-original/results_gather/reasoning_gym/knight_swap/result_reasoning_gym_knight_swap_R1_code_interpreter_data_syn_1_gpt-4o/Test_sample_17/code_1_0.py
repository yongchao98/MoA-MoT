from collections import deque

def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(position, board):
    possible_moves = []
    directions = [(2, 1), (2, -1), (-2, 1), (-2, -1), 
                  (1, 2), (1, -2), (-1, 2), (-1, -2)]
    for dx, dy in directions:
        new_x = chr(ord(position[0]) + dx)
        new_y = str(int(position[1]) + dy)
        new_pos = new_x + new_y
        if new_pos in board and board[new_pos] == '.':
            possible_moves.append(new_pos)
    return possible_moves

def swap_knights_bfs():
    initial_board = {
        'A1': 'B', 'B1': '.', 'C1': 'B', 'D1': '.',
        'A2': 'w', 'B2': '.', 'C2': '.', 'D2': '.',
        'A3': '.', 'B3': '.', 'C3': 'w', 'D3': '.'
    }
    target_board = {
        'A1': 'w', 'C1': 'w',
        'A2': 'B', 'C3': 'B'
    }
    initial_state = ('A2', 'C3', 'A1', 'C1')  # (w1, w2, B1, B2)
    target_state = ('A1', 'C1', 'A2', 'C3')
    
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, moves = queue.popleft()
        w1, w2, B1, B2 = current_state
        
        # Check if target state is reached
        if current_state == target_state:
            print(moves)
            return
        
        # Generate possible moves for white knights
        for w_pos in [w1, w2]:
            possible_moves = get_possible_moves(w_pos, initial_board)
            for move in possible_moves:
                new_state = tuple(sorted([move if x == w_pos else x for x in current_state]))
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, moves + [f"w,{w_pos},{move}"]))
        
        # Generate possible moves for black knights
        for B_pos in [B1, B2]:
            possible_moves = get_possible_moves(B_pos, initial_board)
            for move in possible_moves:
                new_state = tuple(sorted([move if x == B_pos else x for x in current_state]))
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, moves + [f"B,{B_pos},{move}"]))
    
    print("No")

swap_knights_bfs()