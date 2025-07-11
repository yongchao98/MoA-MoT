import math
from functools import lru_cache

# This script solves the described board game by finding the optimal strategy.
# It uses a minimax algorithm with memoization to determine the shortest number of turns
# for Player 1 to force a win against an optimal opponent.

# --- Game Constants ---
BOARD_SIZE = 8
PIECES = ['K1', 'N1', 'R1', 'K2', 'N2', 'R2']
P1_PIECE_INDICES = [0, 1, 2]
P2_PIECE_INDICES = [3, 4, 5]
P1, P2 = 1, 2
P1_WIN, P2_WIN, DRAW = 1, 2, 0

# Initial state: (K1_pos, N1_pos, R1_pos, K2_pos, N2_pos, R2_pos)
initial_state = (0, 1, 2, 7, 6, 5)

def is_check(state, player):
    """Checks if the specified player's King is under attack by the opponent's Rook."""
    king_idx = 0 if player == P1 else 3
    rook_idx = 5 if player == P1 else 2

    king_pos = state[king_idx]
    rook_pos = state[rook_idx]

    if king_pos == -1 or rook_pos == -1:
        return False

    start, end = min(king_pos, rook_pos), max(king_pos, rook_pos)

    # Check for any other pieces blocking the line of sight
    for i, pos in enumerate(state):
        if pos != -1 and i != king_idx and i != rook_idx and start < pos < end:
            return False  # Path is blocked
    return True  # Path is clear, King is in check

def generate_moves(state, player):
    """Generates all legal moves for the given player from the current state."""
    moves = []
    piece_indices = P1_PIECE_INDICES if player == P1 else P2_PIECE_INDICES
    friendly_pos = {state[i] for i in piece_indices if state[i] != -1}
    all_occupied_pos = {pos for pos in state if pos != -1}

    for i in piece_indices:
        pos = state[i]
        if pos == -1: continue

        piece_type = PIECES[i][0]
        potential_dests = []

        if piece_type == 'K':
            potential_dests.extend([pos - 1, pos + 1])
        elif piece_type == 'N':
            potential_dests.extend([pos - 2, pos + 2])
        else:  # Rook
            for new_pos in range(pos + 1, BOARD_SIZE):
                potential_dests.append(new_pos)
                if new_pos in all_occupied_pos: break
            for new_pos in range(pos - 1, -1, -1):
                potential_dests.append(new_pos)
                if new_pos in all_occupied_pos: break

        for new_pos in potential_dests:
            if not (0 <= new_pos < BOARD_SIZE and new_pos not in friendly_pos):
                continue

            new_state_list = list(state)
            new_state_list[i] = new_pos

            if new_pos in all_occupied_pos:
                for j, p_pos in enumerate(state):
                    if p_pos == new_pos and j not in piece_indices:
                        new_state_list[j] = -1
                        break

            new_state = tuple(new_state_list)
            if not is_check(new_state, player):
                moves.append(new_state)
    return moves

def is_better_for_player(new_res, old_res, player):
    """Compares two outcomes to determine the better one for the current player."""
    new_outcome, new_plies = new_res
    old_outcome, old_plies = old_res

    p_wins = P1_WIN if player == P1 else P2_WIN
    p_loses = P2_WIN if player == P1 else P1_WIN

    if new_outcome == p_wins and old_outcome != p_wins: return True
    if new_outcome != p_wins and old_outcome == p_wins: return False
    if new_outcome == p_wins: return new_plies < old_plies

    if new_outcome == DRAW and old_outcome == p_loses: return True
    if new_outcome == p_loses and old_outcome == DRAW: return False

    if new_outcome == old_outcome:
        # If losing (e.g., P1 choosing between two P2_WIN states), choice doesn't matter for win/loss.
        # If P2 choosing between two P1_WIN states, P2 prefers to stall (maximize plies).
        # If drawing, prefer to stall.
        return new_plies > old_plies
    return False

@lru_cache(maxsize=None)
def solve(state, player):
    """
    Recursively solves the game using minimax with memoization (via lru_cache).
    Returns a tuple of (outcome, plies).
    """
    if state[3] == -1: return (P1_WIN, 0)
    if state[0] == -1: return (P2_WIN, 0)

    next_states = generate_moves(state, player)

    if not next_states:
        return (P2_WIN if player == P1 else P1_WIN, 0) if is_check(state, player) else (DRAW, 0)

    best_res = None
    for s_prime in next_states:
        res = solve(s_prime, P2 if player == P1 else P1)
        current_res = (res[0], res[1] + 1)
        if best_res is None or is_better_for_player(current_res, best_res, player):
            best_res = current_res

    return best_res

def find_shortest_forced_win():
    """
    Analyzes the game from the initial state and prints the result.
    """
    outcome, plies = solve(initial_state, P1)
    
    if outcome == P1_WIN:
        # A turn consists of one move by each player.
        # Number of turns is the number of P1's moves, which is ceil(plies / 2).
        turns = math.ceil(plies / 2)
        print(turns)
    elif outcome == P2_WIN:
        print("Player 2 can force a win.")
    else:
        print("The game ends in a draw with optimal play.")

# Run the analysis
find_shortest_forced_win()
