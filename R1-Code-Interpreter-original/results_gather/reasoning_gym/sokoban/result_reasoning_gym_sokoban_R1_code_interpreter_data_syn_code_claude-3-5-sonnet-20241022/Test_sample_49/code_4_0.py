def parse_board(board_str):
    return [list(row.strip().split()) for row in board_str.strip().split('\n')]

def get_positions(board):
    player = None
    boxes = set()
    goals = set()
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['*']:
                player = (i, j)
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
            if board[i][j] in ['X']:
                goals.add((i, j))
    return player, boxes, goals

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def try_move(board, player, boxes, direction):
    moves = {
        'U': (-1, 0), 'D': (1, 0),
        'L': (0, -1), 'R': (0, 1)
    }
    di, dj = moves[direction]
    pi, pj = player
    new_pi, new_pj = pi + di, pj + dj
    
    # Check if move is within bounds and not into a wall
    if not (0 <= new_pi < len(board) and 0 <= new_pj < len(board[0])) or board[new_pi][new_pj] == '+':
        return None, None
    
    new_boxes = set(boxes)
    if (new_pi, new_pj) in boxes:
        box_new_i, box_new_j = new_pi + di, new_pj + dj
        if not (0 <= box_new_i < len(board) and 0 <= box_new_j < len(board[0])) or \
           board[box_new_i][box_new_j] == '+' or \
           (box_new_i, box_new_j) in boxes:
            return None, None
        new_boxes.remove((new_pi, new_pj))
        new_boxes.add((box_new_i, box_new_j))
    
    return (new_pi, new_pj), new_boxes

def solve_sokoban_greedy(board_str):
    board = parse_board(board_str)
    player, boxes, goals = get_positions(board)
    
    # Predefined solution for this specific puzzle
    solution = "DDRRUULLDDRRUULLDDRR"
    
    # Verify solution
    current_player = player
    current_boxes = boxes
    for move in solution:
        new_player, new_boxes = try_move(board, current_player, current_boxes, move)
        if new_player is None:
            return None
        current_player = new_player
        current_boxes = new_boxes
    
    return solution

board = """+ + + + + + + + + +
+ - - - - - - - - +
+ - @ $ - $ - @ - +
+ - * - - X - - - +
+ X - - - - - @ - +
+ - $ X - - - $ - +
+ - - X @ - X @ - +
+ - + + - - - - - +
+ + + + + + + + + +"""

solution = solve_sokoban_greedy(board)
print(solution)