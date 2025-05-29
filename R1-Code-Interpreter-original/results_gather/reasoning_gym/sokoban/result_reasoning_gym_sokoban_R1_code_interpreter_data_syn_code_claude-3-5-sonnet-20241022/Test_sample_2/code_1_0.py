from collections import deque
import heapq

def manhattan_distance(x1, y1, x2, y2):
    return abs(x1 - x2) + abs(y1 - y2)

def get_goals(board):
    goals = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['X', '$']:
                goals.append((i, j))
    return goals

def get_boxes(board):
    boxes = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['@', '$']:
                boxes.append((i, j))
    return boxes

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def is_deadlock(board, box_x, box_y):
    # Check if box is in corner
    if (board[box_x-1][box_y] in ['+'] and board[box_x][box_y-1] in ['+']) or \
       (board[box_x-1][box_y] in ['+'] and board[box_x][box_y+1] in ['+']) or \
       (board[box_x+1][box_y] in ['+'] and board[box_x][box_y-1] in ['+']) or \
       (board[box_x+1][box_y] in ['+'] and board[box_x][box_y+1] in ['+']):
        return True
    return False

def solve_sokoban(board):
    rows, cols = len(board), len(board[0])
    goals = get_goals(board)
    initial_boxes = get_boxes(board)
    initial_player = get_player_pos(board)
    
    # State: (player_pos, boxes, path)
    initial_state = (initial_player, tuple(sorted(initial_boxes)), "")
    visited = {(initial_player, tuple(sorted(initial_boxes)))}
    
    # Priority queue based on heuristic
    pq = [(0, initial_state)]
    moves = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    
    while pq:
        _, (player_pos, boxes, path) = heapq.heappop(pq)
        
        # Check if solved
        if all((box_x, box_y) in goals for box_x, box_y in boxes):
            return path
            
        px, py = player_pos
        for move, dx, dy in moves:
            new_px, new_py = px + dx, py + dy
            
            if board[new_px][new_py] == '+':
                continue
                
            # If moving to box
            if (new_px, new_py) in boxes:
                new_box_x, new_box_y = new_px + dx, new_py + dy
                
                if (board[new_box_x][new_box_y] == '+' or 
                    (new_box_x, new_box_y) in boxes or 
                    is_deadlock(board, new_box_x, new_box_y)):
                    continue
                    
                new_boxes = list(boxes)
                new_boxes.remove((new_px, new_py))
                new_boxes.append((new_box_x, new_box_y))
                new_boxes = tuple(sorted(new_boxes))
                
                if ((new_px, new_py), new_boxes) not in visited:
                    visited.add(((new_px, new_py), new_boxes))
                    # Calculate heuristic (sum of manhattan distances from boxes to goals)
                    h = sum(min(manhattan_distance(bx, by, gx, gy) 
                              for gx, gy in goals)
                          for bx, by in new_boxes)
                    heapq.heappush(pq, (h + len(path) + 1, 
                                      ((new_px, new_py), new_boxes, path + move)))
            else:
                if ((new_px, new_py), boxes) not in visited:
                    visited.add(((new_px, new_py), boxes))
                    h = sum(min(manhattan_distance(bx, by, gx, gy) 
                              for gx, gy in goals)
                          for bx, by in boxes)
                    heapq.heappush(pq, (h + len(path) + 1, 
                                      ((new_px, new_py), boxes, path + move)))
    return None

# Initialize the puzzle
puzzle = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '*', '@', '-', '-', 'X', 'X', 'X', '-', '+'],
    ['+', '@', '$', '@', '-', '@', '-', '@', '-', '+'],
    ['+', 'X', 'X', '-', '-', '-', '-', '@', '-', '+'],
    ['+', '+', '+', '@', '-', 'X', '+', '-', '-', '+'],
    ['+', '+', '+', 'X', '-', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(puzzle)
print(solution if solution else "No solution found")