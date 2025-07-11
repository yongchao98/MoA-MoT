import sys

# Increase recursion limit for deep game trees, although it may not be necessary for this specific problem.
sys.setrecursionlimit(2000)

class GameSolver:
    """
    Solves the 1D chess-like game using a recursive approach with memoization.
    """
    def __init__(self):
        self.BOARD_SIZE = 8
        self.P1 = 1
        self.P2 = 2
        self.DRAW = 0

        # Piece indices in the state tuple for easy access
        self.K1, self.N1, self.R1, self.K2, self.N2, self.R2 = 0, 1, 2, 3, 4, 5

        self.P1_PIECES = [self.K1, self.N1, self.R1]
        self.P2_PIECES = [self.K2, self.N2, self.R2]
        
        self.PIECE_TYPE = {
            self.K1: 'K', self.N1: 'N', self.R1: 'R',
            self.K2: 'K', self.N2: 'N', self.R2: 'R'
        }

        # Initial state: [K1][N1][R1][ ][ ][R2][N2][K2]
        # State tuple: (pos_K1, pos_N1, pos_R1, pos_K2, pos_N2, pos_R2)
        self.initial_state = (0, 1, 2, 7, 6, 5)
        
        # Cache for memoization to store results of (state, player) pairs
        self.cache = {}

    def get_piece_at(self, pos, state):
        """Returns the index of the piece at a given position, or None if empty."""
        for i in range(len(state)):
            if state[i] == pos:
                return i
        return None

    def is_king_in_check(self, king_pos, rook_pos, state):
        """Checks if a King is under attack by a Rook."""
        if king_pos == -1 or rook_pos == -1:
            return False
        
        start = min(king_pos, rook_pos)
        end = max(king_pos, rook_pos)
        
        # Check for any blocking pieces between the King and the Rook
        for i in range(start + 1, end):
            if self.get_piece_at(i, state) is not None:
                return False  # Path is blocked
        
        return True

    def get_moves(self, state, player):
        """Generates all legal moves for a given player from a given state."""
        moves = []
        player_pieces = self.P1_PIECES if player == self.P1 else self.P2_PIECES
        
        occupied_map = {p: i for i, p in enumerate(state) if p != -1}

        for piece_idx in player_pieces:
            pos = state[piece_idx]
            if pos == -1:  # Piece is captured
                continue
            
            piece_type = self.PIECE_TYPE[piece_idx]
            potential_new_pos = []

            if piece_type == 'K':
                potential_new_pos.extend([pos - 1, pos + 1])
            elif piece_type == 'N':
                potential_new_pos.extend([pos - 2, pos + 2])
            elif piece_type == 'R':
                # Move left until blocked or edge
                for p in range(pos - 1, -1, -1):
                    potential_new_pos.append(p)
                    if p in occupied_map:
                        break
                # Move right until blocked or edge
                for p in range(pos + 1, self.BOARD_SIZE):
                    potential_new_pos.append(p)
                    if p in occupied_map:
                        break

            for new_pos in potential_new_pos:
                if not (0 <= new_pos < self.BOARD_SIZE):
                    continue
                
                # Destination cannot be occupied by a friendly piece
                if new_pos in occupied_map and occupied_map[new_pos] in player_pieces:
                    continue
                
                # Generate the next state after the move
                next_state_list = list(state)
                next_state_list[piece_idx] = new_pos
                
                # Handle capture of an opponent's piece
                if new_pos in occupied_map:
                    captured_piece_idx = occupied_map[new_pos]
                    next_state_list[captured_piece_idx] = -1
                
                next_state = tuple(next_state_list)

                # Final legality check: King safety
                my_king_pos = next_state[self.K1] if player == self.P1 else next_state[self.K2]
                opponent_rook_pos = next_state[self.R2] if player == self.P1 else next_state[self.R1]
                
                if not self.is_king_in_check(my_king_pos, opponent_rook_pos, next_state):
                    moves.append(next_state)

        return moves

    def solve(self, state, player):
        """
        Recursively solves the game from the given state and player's turn.
        Returns a tuple (winner, plys_to_end).
        """
        # Return cached result if this state has been seen before
        if (state, player) in self.cache:
            return self.cache[(state, player)]

        # Base case: Win by capturing the opponent's king
        if state[self.K2] == -1: return (self.P1, 0)
        if state[self.K1] == -1: return (self.P2, 0)
        
        possible_moves = self.get_moves(state, player)

        # Base case: Stalemate (no legal moves)
        if not possible_moves:
            return (self.DRAW, 0)

        next_player = self.P2 if player == self.P1 else self.P1
        
        p1_wins, p2_wins, draws = [], [], []

        for move in possible_moves:
            winner, plys = self.solve(move, next_player)
            if winner == self.P1: p1_wins.append(plys)
            elif winner == self.P2: p2_wins.append(plys)
            else: draws.append(plys)
        
        # Apply minimax logic to choose the best outcome
        if player == self.P1:
            if p1_wins: # P1 can win, chooses the fastest path
                result = (self.P1, 1 + min(p1_wins))
            elif draws: # P1 can only draw, stalls for it
                result = (self.DRAW, 1 + max(draws))
            else: # P1 will lose, stalls for as long as possible
                result = (self.P2, 1 + max(p2_wins))
        else: # Player is P2
            if p2_wins: # P2 can win, chooses the fastest path
                result = (self.P2, 1 + min(p2_wins))
            elif draws: # P2 can only draw, stalls for it
                result = (self.DRAW, 1 + max(draws))
            else: # P2 will lose, stalls for as long as possible
                result = (self.P1, 1 + max(p1_wins))

        self.cache[(state, player)] = result
        return result

    def find_solution(self):
        """Runs the solver from the initial state and prints the result."""
        print("Analyzing the game tree to find the optimal strategy...")
        winner, total_plys = self.solve(self.initial_state, self.P1)
        
        print("-" * 30)
        if winner == self.P1:
            # A "turn" for P1 is one of their moves.
            # Turns = ceil(plys / 2). (plys + 1) // 2 is integer ceil.
            turns = (total_plys + 1) // 2
            print(f"Outcome: Player 1 can force a win.")
            print(f"Total plys (half-turns) required: {total_plys}")
            print(f"This means the win is achieved on Player 1's turn number: {turns}")
        elif winner == self.P2:
            turns = (total_plys + 1) // 2
            print(f"Outcome: Player 2 can force a win in {turns} of their turns.")
        else:
            print("Outcome: The game is a forced draw.")

if __name__ == '__main__':
    solver = GameSolver()
    solver.find_solution()
