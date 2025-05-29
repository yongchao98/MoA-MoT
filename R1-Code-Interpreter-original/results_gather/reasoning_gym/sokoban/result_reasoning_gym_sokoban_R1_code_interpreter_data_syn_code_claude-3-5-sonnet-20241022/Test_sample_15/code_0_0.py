from collections import deque
import copy

def parse_board(board):
    rows = board.strip().split('\n')
    height = len(rows)
    width = len(rows[0])
    player = None
    boxes = set()
    goals = set()
    walls = set()
    
    for y, row in enumerate(rows):
        for x, cell in enumerate(row.strip()):
            pos = (x, y)
            if cell in ['*', '%']:  # Player
                player = pos
            if cell in ['@', '$']:  # Box
                boxes.add(pos)
            if cell in ['X', '$', '%']:  # Goal
                goals.add(pos)
            if cell == '+':  # Wall
                walls.add(pos)
                
    return player, boxes, goals, walls, width, height

def get_moves(player, boxes, walls, width, height):
    moves = []
    directions = {'U': (0, -1), 'D': (0, 1), 'L': (-1, 0), 'R': (1, 0)}
    
    for move, (dx, dy) in directions.items():
        new_x, new_y = player[0] + dx, player[1] + dy
        if 0 <= new_x < width and 0 <= new_y < height:
            new_pos = (new_x, new_y)
            if new_pos not in walls:
                if new_pos in boxes:
                    box_x, box_y = new_x + dx, new_y + dy
                    if (box_x, box_y) not in walls and (box_x, box_y) not in boxes:
                        moves.append((move, new_pos, (box_x, box_y)))
                else:
                    moves.append((move, new_pos, None))
    return moves

def solve_sokoban(board):
    player, boxes, goals, walls, width, height = parse_board(board)
    queue = deque([(player, frozenset(boxes), "")])
    visited = set()
    
    while queue:
        player, boxes, path = queue.popleft()
        
        if boxes == goals:
            return path
            
        state = (player, boxes)
        if state in visited:
            continue
        visited.add(state)
        
        for move, new_player, new_box in get_moves(player, boxes, walls, width, height):
            new_boxes = set(boxes)
            if new_box:
                new_boxes.remove(new_player)
                new_boxes.add(new_box)
            new_state = (new_player, frozenset(new_boxes))
            if new_state not in visited:
                queue.append((new_player, frozenset(new_boxes), path + move))
    
    return None

# Test the puzzle
puzzle = """+ + + + + + + +
+ X + $ X @ * +
+ - - - - - @ +
+ @ - $ - - - +
+ - - - - - X +
+ + + + + + + +"""

solution = solve_sokoban(puzzle)
print(solution)