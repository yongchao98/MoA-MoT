from collections import deque, defaultdict
import heapq

def parse_board(board_str):
    lines = board_str.strip().split('\n')
    board = [list(line.strip()) for line in lines]
    player = None
    boxes = set()
    goals = set()
    
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '%']:
                player = (i, j)
                if board[i][j] == '%':
                    goals.add((i, j))
            elif board[i][j] == '$':
                boxes.add((i, j))
                goals.add((i, j))
            elif board[i][j] == 'X':
                goals.add((i, j))
    
    return board, player, boxes, goals

def is_valid_move(board, pos):
    return 0 <= pos[0] < len(board) and 0 <= pos[1] < len(board[0]) and board[pos[0]][pos[1]] != '+'

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def get_min_matching(boxes, goals):
    total = 0
    unmatched_boxes = list(boxes)
    unmatched_goals = list(goals)
    
    while unmatched_boxes:
        min_dist = float('inf')
        min_box = None
        min_goal = None
        
        for box in unmatched_boxes:
            for goal in unmatched_goals:
                dist = manhattan_distance(box, goal)
                if dist < min_dist:
                    min_dist = dist
                    min_box = box
                    min_goal = goal
        
        total += min_dist
        unmatched_boxes.remove(min_box)
        unmatched_goals.remove(min_goal)
    
    return total

def get_neighbors(pos):
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    return [(pos[0] + dy, pos[1] + dx, move) for dy, dx, move in directions]

def solve_sokoban(board_str):
    board, player, boxes, goals = parse_board(board_str)
    initial_state = (player, frozenset(boxes))
    visited = set()
    pq = [(0, 0, '', initial_state)]
    g_score = defaultdict(lambda: float('inf'))
    g_score[initial_state] = 0
    
    while pq:
        _, cost, path, state = heapq.heappop(pq)
        player, boxes = state
        
        if boxes == goals:
            return path
            
        if state in visited:
            continue
        visited.add(state)
        
        for next_pos, box_pos, move in get_neighbors(player):
            if not is_valid_move(board, next_pos):
                continue
                
            new_boxes = set(boxes)
            if next_pos in boxes:
                box_next_pos = (next_pos[0] + (next_pos[0] - player[0]),
                              next_pos[1] + (next_pos[1] - player[1]))
                if not is_valid_move(board, box_next_pos) or box_next_pos in boxes:
                    continue
                new_boxes.remove(next_pos)
                new_boxes.add(box_next_pos)
            
            new_state = (next_pos, frozenset(new_boxes))
            new_cost = g_score[state] + 1
            
            if new_cost < g_score[new_state]:
                g_score[new_state] = new_cost
                h_score = get_min_matching(new_boxes, goals)
                f_score = new_cost + h_score
                heapq.heappush(pq, (f_score, new_cost, path + move, new_state))
    
    return None

# Define the puzzle
puzzle = """+ + + + + + +
+ $ + - - X +
+ - - - - - +
+ - - - @ - +
+ - X - @ - +
+ - - - X - +
+ - - @ X X +
+ - - - @ @ +
+ - - - @ % +
+ + + + + + +"""

# Solve the puzzle
solution = solve_sokoban(puzzle)
print(solution)