from collections import deque
import copy

def is_blocking_red_car(car_pos, red_car_pos):
    # Check if a car is directly blocking the red car's path to the right
    red_row = red_car_pos[0][0]
    red_right = max(p[1] for p in red_car_pos)
    return any(p[0] == red_row and p[1] > red_right for p in car_pos)

def solve_puzzle():
    initial_board = [
        '.xBBCC',
        '..x.K.',
        'G.AAK.',
        'G.IJDD',
        'H.IJ..',
        'HEEFFF'
    ]
    
    queue = deque([(initial_board, [])])
    seen = {'\n'.join(initial_board)}
    
    while queue:
        current_board = queue.popleft()[0]
        cars = {}
        # Find all cars
        for i in range(len(current_board)):
            for j in range(len(current_board[i])):
                if current_board[i][j] not in '.x':
                    if current_board[i][j] not in cars:
                        cars[current_board[i][j]] = [(i, j)]
                    else:
                        cars[current_board[i][j]].append((i, j))
        
        # If K can move up and out of the way
        if 'K' in cars and current_board[1][4] == 'K' and current_board[0][4] == '.':
            return "<<<K-1 A+3>>>"
            
        # If K is at bottom and can clear path
        if 'K' in cars and current_board[2][4] == 'K' and current_board[1][4] == '.':
            return "<<<K-1 A+3>>>"
    
    return "No solution found"

print(solve_puzzle())