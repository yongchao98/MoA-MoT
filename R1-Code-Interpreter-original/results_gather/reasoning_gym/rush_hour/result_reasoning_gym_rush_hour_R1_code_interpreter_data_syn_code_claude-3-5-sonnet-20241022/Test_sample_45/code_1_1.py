def solve_rush_hour(board):
    initial_board = parse_board(board)
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        vehicles = get_vehicle_info(current_board)
        
        # Check if red car (A) can reach exit
        a_coords = vehicles['A']
        if max(c[1] for c in a_coords) == len(board[0]) - 1:
            return moves
        
        possible_moves = get_moves(current_board, vehicles)
        
        for vehicle, direction in possible_moves:
            new_board, new_coords = apply_move(current_board, vehicles, vehicle, direction)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                move_str = f"{vehicle}{'+' if direction > 0 else ''}{direction}"
                queue.append((new_board, moves + [move_str]))
    
    return None

# Initial board
board = [
    "BBH.CC",
    ".GH..K",
    ".GAAJK",
    "DD.IJL",
    "...IxL",
    ".EEFFL"
]

solution = solve_rush_hour(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")