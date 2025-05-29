def is_deadlock(board, box_pos):
    row, col = box_pos
    
    # Corner deadlock check
    if (board[row-1][col] == '+' and board[row][col-1] == '+') or \
       (board[row-1][col] == '+' and board[row][col+1] == '+') or \
       (board[row+1][col] == '+' and board[row][col-1] == '+') or \
       (board[row+1][col] == '+' and board[row][col+1] == '+'):
        return True
    return False

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def solve_sokoban(board_str):
    # Parse board
    board = [list(row.strip().split()) for row in board_str.strip().split('\n')]
    rows, cols = len(board), len(board[0])
    
    # Find initial positions
    player = None
    boxes = []
    goals = []
    for i in range(rows):
        for j in range(cols):
            if board[i][j] in ['*', '%']:
                player = (i, j)
            elif board[i][j] == '$':
                boxes.append((i, j))
            elif board[i][j] == 'X':
                goals.append((i, j))
            elif board[i][j] == '@':
                boxes.append((i, j))
                goals.append((i, j))

    def get_moves(pos, boxes_set):
        moves = []
        for dx, dy, move in [(-1,0,'U'), (1,0,'D'), (0,-1,'L'), (0,1,'R')]:
            new_pos = (pos[0] + dx, pos[1] + dy)
            
            # Check if move is valid
            if board[new_pos[0]][new_pos[1]] == '+':
                continue
                
            if new_pos in boxes_set:
                push_pos = (new_pos[0] + dx, new_pos[1] + dy)
                if board[push_pos[0]][push_pos[1]] != '+' and push_pos not in boxes_set:
                    if not is_deadlock(board, push_pos):
                        moves.append((move, new_pos, push_pos))
            else:
                moves.append((move, new_pos, None))
        return moves

    def solve_dfs(player, boxes, path, depth=0, max_depth=50):
        if depth >= max_depth:
            return None
            
        # Check if won
        if all(box in goals for box in boxes):
            return path
            
        boxes_set = set(boxes)
        for move, new_player, push_pos in get_moves(player, boxes_set):
            if push_pos:
                new_boxes = tuple(push_pos if b == new_player else b for b in boxes)
            else:
                new_boxes = boxes
                
            result = solve_dfs(new_player, new_boxes, path + move, depth + 1)
            if result:
                return result
        return None

    # Try increasing depth limits
    for max_depth in [20, 30, 40, 50]:
        solution = solve_dfs(player, tuple(boxes), "", max_depth=max_depth)
        if solution:
            return solution
    return None

puzzle = """+ + + + + + +
+ - - - - - +
+ $ @ X - - +
+ - - - @ - +
+ - - - $ - +
+ X - - - - +
+ + @ X @ - +
+ + X - - * +
+ + + + + + +"""

print(solve_sokoban(puzzle))