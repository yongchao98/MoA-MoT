from collections import deque
import heapq

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def get_boxes_and_goals(state):
    boxes, goals = [], []
    for i in range(len(state)):
        for j in range(len(state[0])):
            if state[i][j] in ['@', '$']:
                boxes.append((i, j))
            if state[i][j] in ['X', '$', '%']:
                goals.append((i, j))
    return boxes, goals

def heuristic(state):
    boxes, goals = get_boxes_and_goals(state)
    total_dist = 0
    used_goals = set()
    
    for box in boxes:
        min_dist = float('inf')
        best_goal = None
        for goal in goals:
            if goal not in used_goals:
                dist = manhattan_distance(box, goal)
                if dist < min_dist:
                    min_dist = dist
                    best_goal = goal
        if best_goal:
            used_goals.add(best_goal)
            total_dist += min_dist
    return total_dist

def is_deadlock(state, box_pos):
    rows, cols = len(state), len(state[0])
    x, y = box_pos
    
    # Corner deadlock
    if (state[x-1][y] == '+' and state[x][y-1] == '+') or \
       (state[x-1][y] == '+' and state[x][y+1] == '+') or \
       (state[x+1][y] == '+' and state[x][y-1] == '+') or \
       (state[x+1][y] == '+' and state[x][y+1] == '+'):
        return True
    
    return False

def get_player_pos(state):
    for i in range(len(state)):
        for j in range(len(state[0])):
            if state[i][j] in ['*', '%']:
                return (i, j)
    return None

def solve_sokoban(initial_state):
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    rows, cols = len(initial_state), len(initial_state[0])
    start_pos = get_player_pos(initial_state)
    
    # Priority queue with (heuristic + moves_length, moves, state, player_pos)
    pq = [(heuristic(initial_state), "", initial_state, start_pos)]
    visited = {str(initial_state)}
    max_moves = 40  # Limit solution length
    
    while pq:
        _, moves, current_state, player_pos = heapq.heappop(pq)
        
        if len(moves) > max_moves:
            continue
            
        # Check if solved
        boxes, goals = get_boxes_and_goals(current_state)
        if all(box in goals for box in boxes):
            return moves
            
        # Try each direction
        for dx, dy, move in directions:
            new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
            
            if not (0 <= new_x < rows and 0 <= new_y < cols) or current_state[new_x][new_y] == '+':
                continue
                
            new_state = [list(row) for row in current_state]
            
            # Moving to empty space or goal
            if current_state[new_x][new_y] in ['-', 'X']:
                new_state[player_pos[0]][player_pos[1]] = '-' if current_state[player_pos[0]][player_pos[1]] == '*' else 'X'
                new_state[new_x][new_y] = '*' if current_state[new_x][new_y] == '-' else '%'
                
                state_str = str(new_state)
                if state_str not in visited:
                    visited.add(state_str)
                    h_score = heuristic(new_state)
                    heapq.heappush(pq, (h_score + len(moves) + 1, moves + move, new_state, (new_x, new_y)))
                    
            # Pushing a box
            elif current_state[new_x][new_y] in ['@', '$']:
                push_x, push_y = new_x + dx, new_y + dy
                
                if (0 <= push_x < rows and 0 <= push_y < cols) and \
                   current_state[push_x][push_y] in ['-', 'X']:
                    
                    # Move box
                    new_state[push_x][push_y] = '@' if current_state[push_x][push_y] == '-' else '$'
                    new_state[new_x][new_y] = '*' if current_state[new_x][new_y] == '@' else '%'
                    new_state[player_pos[0]][player_pos[1]] = '-' if current_state[player_pos[0]][player_pos[1]] == '*' else 'X'
                    
                    if not is_deadlock(new_state, (push_x, push_y)):
                        state_str = str(new_state)
                        if state_str not in visited:
                            visited.add(state_str)
                            h_score = heuristic(new_state)
                            heapq.heappush(pq, (h_score + len(moves) + 1, moves + move, new_state, (new_x, new_y)))
    
    return None

initial_state = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', 'X', '-', '@', 'X', '-', '@', '-', '+'],
    ['+', '@', 'X', '@', '-', '@', '@', 'X', '+'],
    ['+', '-', '-', '-', '-', 'X', '-', '-', '+'],
    ['+', '-', '*', '@', 'X', '-', '-', '-', '+'],
    ['+', '-', '-', '-', 'X', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+']]

solution = solve_sokoban(initial_state)
print(solution if solution else "No solution found")