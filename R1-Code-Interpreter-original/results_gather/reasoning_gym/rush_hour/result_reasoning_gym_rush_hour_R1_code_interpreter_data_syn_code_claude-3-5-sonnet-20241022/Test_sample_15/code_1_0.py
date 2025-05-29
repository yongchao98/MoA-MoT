from heapq import heappush, heappop
import copy

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def heuristic(board, vehicles):
    # Distance of AA from exit plus number of blocking vehicles
    aa_pos = vehicles['A']
    exit_pos = (2, 5)  # The exit position
    blocking_vehicles = 0
    
    # Count vehicles blocking AA's path to exit
    right_most_aa = max(p[1] for p in aa_pos)
    for col in range(right_most_aa + 1, 6):
        if board[2][col] not in '.x':
            blocking_vehicles += 1
            
    return manhattan_distance(aa_pos[0], (2, 5)) + blocking_vehicles * 2

def solve_puzzle(board_str):
    board = [list(row) for row in board_str.split('\n') if row]
    vehicles = {}
    
    # Initialize vehicles positions
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in '.x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = [(i, j)]
                else:
                    vehicles[board[i][j]].append((i, j))
    
    # Priority queue for A* search: (priority, moves_count, board_state, moves, vehicles)
    queue = [(heuristic(board, vehicles), 0, board, [], vehicles)]
    seen = {'\n'.join(''.join(row) for row in board)}
    
    while queue:
        _, moves_count, current_board, moves, current_vehicles = heappop(queue)
        
        # Check if solved
        aa_pos = current_vehicles['A']
        if max(p[1] for p in aa_pos) == 5:
            return moves
        
        # Generate possible moves
        for vehicle, positions in current_vehicles.items():
            is_horizontal = positions[0][0] == positions[1][0]
            
            if is_horizontal:
                row = positions[0][0]
                min_col = min(p[1] for p in positions)
                max_col = max(p[1] for p in positions)
                
                # Try moving left
                if min_col > 0 and current_board[row][min_col-1] == '.':
                    new_board = [row[:] for row in current_board]
                    for pos in positions:
                        new_board[pos[0]][pos[1]] = '.'
                        new_board[pos[0]][pos[1]-1] = vehicle
                    
                    board_str = '\n'.join(''.join(row) for row in new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        new_vehicles = {k: [(p[0], p[1]-1) if k==vehicle else p for p in v] 
                                     for k,v in current_vehicles.items()}
                        priority = moves_count + 1 + heuristic(new_board, new_vehicles)
                        heappush(queue, (priority, moves_count + 1, new_board, 
                                       moves + [(vehicle, -1)], new_vehicles))
                
                # Try moving right
                if max_col < 5 and current_board[row][max_col+1] == '.':
                    new_board = [row[:] for row in current_board]
                    for pos in positions:
                        new_board[pos[0]][pos[1]] = '.'
                        new_board[pos[0]][pos[1]+1] = vehicle
                    
                    board_str = '\n'.join(''.join(row) for row in new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        new_vehicles = {k: [(p[0], p[1]+1) if k==vehicle else p for p in v] 
                                     for k,v in current_vehicles.items()}
                        priority = moves_count + 1 + heuristic(new_board, new_vehicles)
                        heappush(queue, (priority, moves_count + 1, new_board, 
                                       moves + [(vehicle, 1)], new_vehicles))
            else:
                col = positions[0][1]
                min_row = min(p[0] for p in positions)
                max_row = max(p[0] for p in positions)
                
                # Try moving up
                if min_row > 0 and current_board[min_row-1][col] == '.':
                    new_board = [row[:] for row in current_board]
                    for pos in positions:
                        new_board[pos[0]][pos[1]] = '.'
                        new_board[pos[0]-1][pos[1]] = vehicle
                    
                    board_str = '\n'.join(''.join(row) for row in new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        new_vehicles = {k: [(p[0]-1, p[1]) if k==vehicle else p for p in v] 
                                     for k,v in current_vehicles.items()}
                        priority = moves_count + 1 + heuristic(new_board, new_vehicles)
                        heappush(queue, (priority, moves_count + 1, new_board, 
                                       moves + [(vehicle, -1)], new_vehicles))
                
                # Try moving down
                if max_row < 5 and current_board[max_row+1][col] == '.':
                    new_board = [row[:] for row in current_board]
                    for pos in positions:
                        new_board[pos[0]][pos[1]] = '.'
                        new_board[pos[0]+1][pos[1]] = vehicle
                    
                    board_str = '\n'.join(''.join(row) for row in new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        new_vehicles = {k: [(p[0]+1, p[1]) if k==vehicle else p for p in v] 
                                     for k,v in current_vehicles.items()}
                        priority = moves_count + 1 + heuristic(new_board, new_vehicles)
                        heappush(queue, (priority, moves_count + 1, new_board, 
                                       moves + [(vehicle, 1)], new_vehicles))

board = """.IBBBx
GI.JCC
GAAJK.
G.DDKL
H.xEEL
HFFF.."""

solution = solve_puzzle(board)
if solution:
    formatted_moves = []
    for vehicle, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_moves.append(f"{vehicle}{sign}1")
    print("<<<" + " ".join(formatted_moves) + ">>>")
else:
    print("<<<No solution found>>>")