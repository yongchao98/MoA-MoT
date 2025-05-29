import heapq

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def get_positions(board, chars):
    positions = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in chars:
                positions.append((i, j))
    return positions

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_sokoban(board):
    rows, cols = len(board), len(board[0])
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    # Get initial positions
    player = next((i, j) for i in range(rows) for j in range(cols) if board[i][j] in ['*', '%'])
    boxes = [(i, j) for i in range(rows) for j in range(cols) if board[i][j] in ['@', '$']]
    goals = [(i, j) for i in range(rows) for j in range(cols) if board[i][j] in ['X', '$', '%']]
    
    def heuristic(player_pos, box_positions):
        # Minimum distance from boxes to goals
        total = 0
        unused_goals = goals.copy()
        for box in sorted(box_positions):
            if unused_goals:
                min_dist = min(manhattan_distance(box, goal) for goal in unused_goals)
                total += min_dist
                # Remove the closest goal
                closest_goal = min(unused_goals, key=lambda g: manhattan_distance(box, g))
                unused_goals.remove(closest_goal)
        return total

    initial_state = (player, tuple(sorted(boxes)))
    visited = {initial_state}
    pq = [(0, 0, '', initial_state)]
    max_steps = 100  # Limit search depth
    
    while pq and max_steps > 0:
        max_steps -= 1
        _, steps, path, (curr_player, curr_boxes) = heapq.heappop(pq)
        
        # Check if solved
        if all(box in goals for box in curr_boxes):
            return path
        
        # Try each direction
        for move, (dx, dy) in directions.items():
            new_player = (curr_player[0] + dx, curr_player[1] + dy)
            
            # Check if player move is valid
            if board[new_player[0]][new_player[1]] == '+':
                continue
                
            new_boxes = list(curr_boxes)
            # If pushing a box
            if new_player in curr_boxes:
                box_idx = new_boxes.index(new_player)
                new_box_pos = (new_player[0] + dx, new_player[1] + dy)
                
                # Check if box move is valid
                if (board[new_box_pos[0]][new_box_pos[1]] == '+' or 
                    new_box_pos in new_boxes):
                    continue
                    
                new_boxes[box_idx] = new_box_pos
            
            new_boxes = tuple(sorted(new_boxes))
            new_state = (new_player, new_boxes)
            
            if new_state not in visited:
                visited.add(new_state)
                h = heuristic(new_player, new_boxes)
                heapq.heappush(pq, (steps + 1 + h, steps + 1, path + move, new_state))
    
    return None

# Initial board
board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '-', '+', '+', '$', '+', '+', '+'],
    ['+', 'X', '-', '@', '-', '-', 'X', '$', '+'],
    ['+', '@', '-', '+', '$', 'X', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '+', 'X', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '-', 'X', '-', '$', '-', '-', '-', '+'],
    ['+', '-', '@', '-', '@', '-', '-', '@', '+'],
    ['+', '-', '-', '-', '*', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(board)
print(solution if solution else "No solution found")