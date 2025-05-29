from collections import deque

def is_valid_move(start, end, occupied):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx, dy) in [(2, 1), (1, 2)] and end not in occupied

def swap_knights():
    # Initial positions
    white_knights = ['B1', 'D1']
    black_knights = ['C3', 'B2']
    empty_squares = ['A3', 'A2', 'C1']
    
    # Target positions
    target_white = ['C3', 'B2']
    target_black = ['B1', 'D1']
    
    # Queue for BFS
    queue = deque()
    queue.append((white_knights, black_knights, empty_squares, [], 'w'))
    
    while queue:
        w_knights, b_knights, empty, moves, turn = queue.popleft()
        
        if set(w_knights) == set(target_white) and set(b_knights) == set(target_black):
            return moves
        
        if turn == 'w':
            for i, wk in enumerate(w_knights):
                for es in empty:
                    if is_valid_move(wk, es, b_knights + empty):
                        new_w_knights = w_knights[:i] + [es] + w_knights[i+1:]
                        new_empty = empty[:]
                        new_empty.remove(es)
                        new_empty.append(wk)
                        new_moves = moves + [f"w,{wk},{es}"]
                        queue.append((new_w_knights, b_knights, new_empty, new_moves, 'B'))
        else:
            for i, bk in enumerate(b_knights):
                for es in empty:
                    if is_valid_move(bk, es, w_knights + empty):
                        new_b_knights = b_knights[:i] + [es] + b_knights[i+1:]
                        new_empty = empty[:]
                        new_empty.remove(es)
                        new_empty.append(bk)
                        new_moves = moves + [f"B,{bk},{es}"]
                        queue.append((w_knights, new_b_knights, new_empty, new_moves, 'w'))
    
    return "No"

# Execute the function and print the result
print(swap_knights())