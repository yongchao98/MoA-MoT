from collections import deque
import heapq
from copy import deepcopy

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def get_positions(board, chars):
    positions = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in chars:
                positions.append((i, j))
    return positions

def is_valid_move(board, player_pos, direction):
    rows, cols = len(board), len(board[0])
    new_r, new_c = player_pos[0] + direction[0], player_pos[1] + direction[1]
    
    if not (0 <= new_r < rows and 0 <= new_c < cols):
        return False
    
    if board[new_r][new_c] == '+':
        return False
        
    if board[new_r][new_c] in '@$':
        box_r, box_c = new_r + direction[0], new_c + direction[1]
        if not (0 <= box_r < rows and 0 <= box_c < cols):
            return False
        if board[box_r][box_c] in '+@$':
            return False
    return True

def make_move(board, player_pos, direction):
    new_board = [list(row) for row in board]
    r, c = player_pos
    new_r, new_c = r + direction[0], c + direction[1]
    
    # Update player position
    if board[r][c] == '*':
        new_board[r][c] = 'X'
    else:
        new_board[r][c] = '-'
        
    if board[new_r][new_c] in '@$':
        box_r, box_c = new_r + direction[0], new_c + direction[1]
        if board[new_r][new_c] == '$':
            new_board[new_r][new_c] = 'X'
        else:
            new_board[new_r][new_c] = '-'
            
        if board[box_r][box_c] == 'X':
            new_board[box_r][box_c] = '$'
        else:
            new_board[box_r][box_c] = '@'
    
    if board[new_r][new_c] == 'X':
        new_board[new_r][new_c] = '*'
    else:
        new_board[new_r][new_c] = '*'
        
    return [''.join(row) for row in new_board], (new_r, new_c)

def heuristic(board):
    boxes = get_positions(board, '@$')
    goals = get_positions(board, 'X*$')
    
    total_dist = 0
    for box in boxes:
        min_dist = float('inf')
        for goal in goals:
            dist = manhattan_distance(box, goal)
            min_dist = min(min_dist, dist)
        total_dist += min_dist
    return total_dist

def solve_sokoban(initial_board):
    directions = {
        'U': (-1, 0), 'D': (1, 0),
        'L': (0, -1), 'R': (0, 1)
    }
    
    player_pos = get_positions(initial_board, '*')[0]
    initial_state = (tuple(initial_board), player_pos)
    visited = set([initial_state])
    pq = [(0, '', initial_state)]
    
    while pq:
        _, path, state = heapq.heappop(pq)
        board, pos = state
        
        if all(c not in ''.join(board) for c in '@'):
            return path
            
        for move, dir_vec in directions.items():
            if is_valid_move(board, pos, dir_vec):
                new_board, new_pos = make_move(board, pos, dir_vec)
                new_state = (tuple(new_board), new_pos)
                
                if new_state not in visited:
                    visited.add(new_state)
                    priority = len(path) + 1 + heuristic(new_board)
                    heapq.heappush(pq, (priority, path + move, new_state))
    
    return None

# Initial board
initial_board = [
    "+++++++" ,
    "+-X--X+",
    "+-@@--+",
    "+X@-X-+",
    "+---+++",
    "+---$+",
    "+-@-+++",
    "+---$+",
    "+-X-@*+",
    "+++++++"
]

solution = solve_sokoban(initial_board)
print(solution)