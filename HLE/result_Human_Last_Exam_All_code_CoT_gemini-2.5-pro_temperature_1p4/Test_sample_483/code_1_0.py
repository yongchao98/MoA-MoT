import sys

# Increase recursion limit for deep searches, just in case.
sys.setrecursionlimit(2000)

# Global memoization table to store results of subproblems
memo = {}

def get_board_from_pos(positions):
    """Converts a dictionary of piece positions to a list representation of the board."""
    board = [' '] * 8
    for piece, pos in positions.items():
        if pos is not None:
            board[pos] = piece
    return board

def is_king_in_check(board, king_player):
    """Checks if the specified player's king is under attack by the opponent's rook."""
    if king_player == 1:
        king_piece, rook_piece = 'K1', 'R2'
    else:
        king_piece, rook_piece = 'K2', 'R1'

    try:
        king_pos = board.index(king_piece)
    except ValueError:
        return False  # King is captured, game over

    try:
        rook_pos = board.index(rook_piece)
    except ValueError:
        return False  # Opponent's rook is captured

    # Check for any pieces between the king and the rook
    start, end = sorted((king_pos, rook_pos))
    for i in range(start + 1, end):
        if board[i] != ' ':
            return False  # Path is blocked
    return True

def get_piece_info(piece_char):
    """Extracts player number and piece type from its string representation."""
    player = int(piece_char[1])
    piece_type = piece_char[0]
    return player, piece_type

def generate_moves(positions, player):
    """Generates all legal moves for a given player from a given position."""
    legal_next_positions = []
    player_pieces = [p for p in positions if p[1] == str(player) and positions.get(p) is not None]
    board = get_board_from_pos(positions)

    for piece in player_pieces:
        current_pos = positions[piece]
        piece_type = piece[0]

        potential_dests = []
        if piece_type == 'K':
            potential_dests.extend([current_pos - 1, current_pos + 1])
        elif piece_type == 'N':
            potential_dests.extend([current_pos - 2, current_pos + 2])
        elif piece_type == 'R':
            for i in range(current_pos - 1, -1, -1):
                potential_dests.append(i)
                if board[i] != ' ': break
            for i in range(current_pos + 1, 8):
                potential_dests.append(i)
                if board[i] != ' ': break

        for dest in potential_dests:
            if not (0 <= dest <= 7): continue
            target_content = board[dest]
            if target_content != ' ' and get_piece_info(target_content)[0] == player: continue

            next_positions = positions.copy()
            next_positions[piece] = dest
            if target_content != ' ':
                next_positions[target_content] = None

            next_board = get_board_from_pos(next_positions)
            if not is_king_in_check(next_board, player):
                legal_next_positions.append(next_positions)
                
    return legal_next_positions

def find_win_in_x_moves(positions, player):
    """
    Recursively solves the game using minimax, returning the minimum number of moves for P1 to win.
    Returns a tuple: (number of moves to win, best next state).
    """
    pos_tuple = tuple(sorted(item for item in positions.items() if item[1] is not None))
    state_key = (pos_tuple, player)

    if state_key in memo: return memo[state_key]

    if positions.get('K2') is None: return (0, None)
    if positions.get('K1') is None: return (float('inf'), None)

    legal_moves = generate_moves(positions, player)
    
    if not legal_moves:
        board = get_board_from_pos(positions)
        if is_king_in_check(board, player):
            return (0, None) if player == 2 else (float('inf'), None)
        else:
            return (float('inf'), None)

    if player == 1:
        best_outcome = (float('inf'), None)
        for next_pos in legal_moves:
            moves, _ = find_win_in_x_moves(next_pos, 2)
            if moves < best_outcome[0]:
                best_outcome = (moves, next_pos)
        result = (1 + best_outcome[0], best_outcome[1]) if best_outcome[0] != float('inf') else (float('inf'), None)
    else:
        best_outcome = (-1, None)
        for next_pos in legal_moves:
            moves, _ = find_win_in_x_moves(next_pos, 1)
            if moves > best_outcome[0]:
                best_outcome = (moves, next_pos)
        result = (1 + best_outcome[0], best_outcome[1]) if best_outcome[0] != float('inf') else (float('inf'), best_outcome[1])

    memo[state_key] = result
    return result

def describe_move(p_from, p_to, player):
    """Generates a human-readable string describing the move between two position states."""
    moved_piece, old_pos, new_pos, captured_piece = None, -1, -1, None
    
    moved_piece = next((p for p, pos in p_from.items() if pos is not None and p_to.get(p) != pos), None)
    if moved_piece:
        old_pos = p_from[moved_piece]
        new_pos = p_to[moved_piece]
        captured_piece = next((p for p, pos in p_from.items() if pos == new_pos and p != moved_piece), None)

    player_name = f"Player {player}"
    move_str = f"{player_name}: {moved_piece} from {old_pos} to {new_pos}"
    if captured_piece:
        move_str += f" (captures {captured_piece})"
    
    return move_str

def solve_game():
    """Main function to solve the game and print the results."""
    initial_positions = {'K1': 0, 'N1': 1, 'R1': 2, 'R2': 5, 'N2': 6, 'K2': 7}

    total_moves, next_pos = find_win_in_x_moves(initial_positions, 1)

    if total_moves == float('inf'):
        print("Player 1 cannot force a win.")
        return

    print(f"Player 1 can force a win in {total_moves} turns.")
    print("The winning sequence of moves is:")
    
    current_pos = initial_positions
    player = 1
    path_str = []
    
    for i in range(total_moves):
        _, best_next_pos = find_win_in_x_moves(current_pos, player)
        move_desc = describe_move(current_pos, best_next_pos, player)
        print(f"Turn {i+1}. {move_desc}")
        path_str.append(str(1)) # Representing one move in the equation
        current_pos = best_next_pos
        player = 2 if player == 1 else 1

    print(f"\nThis corresponds to the equation: {' + '.join(path_str)} = {total_moves}")


if __name__ == '__main__':
    solve_game()
    print("<<<1>>>")
