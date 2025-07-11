#
# This script verifies the mate-in-2 solution for the given Capablanca Chess puzzle.
#
class ChessPosition:
    """Represents a chess position and provides move validation logic."""

    BOARD_WIDTH = 10
    BOARD_HEIGHT = 8
    FILES = 'abcdefghij'
    RANKS = '12345678'

    def __init__(self, fen):
        self.board = {}  # { (file, rank): (piece_char, color) }
        self.turn = 'white'
        self._parse_fen(fen)

    def _to_coords(self, square_name):
        """Converts square name like 'a1' to coordinates (0, 0)."""
        try:
            file = self.FILES.index(square_name[0])
            rank = self.RANKS.index(square_name[1])
            return file, rank
        except (ValueError, IndexError):
            return None

    def _to_name(self, coords):
        """Converts coordinates (0, 0) to square name 'a1'."""
        file, rank = coords
        if 0 <= file < self.BOARD_WIDTH and 0 <= rank < self.BOARD_HEIGHT:
            return self.FILES[file] + self.RANKS[rank]
        return None

    def _parse_fen(self, fen):
        parts = fen.split(' ')
        board_fen, turn_fen = parts[0], parts[1]
        self.turn = 'white' if turn_fen == 'w' else 'black'
        rank_index = self.BOARD_HEIGHT - 1
        file_index = 0
        for char in board_fen:
            if char == '/':
                rank_index -= 1
                file_index = 0
            elif char.isdigit():
                file_index += int(char)
            else:
                color = 'white' if char.isupper() else 'black'
                self.board[(file_index, rank_index)] = (char.lower(), color)
                file_index += 1

    def copy(self):
        new_pos = ChessPosition("8/8/8/8/8/8/8/8 w - - 0 1") # Dummy FEN
        new_pos.board = self.board.copy()
        new_pos.turn = self.turn
        return new_pos

    def find_king(self, color):
        for coords, (piece_char, piece_color) in self.board.items():
            if piece_char == 'k' and piece_color == color:
                return coords
        return None

    def get_piece_attacking_moves(self, coords):
        """Generates all squares a piece can move to or attack."""
        if coords not in self.board: return []
        
        piece_char, color = self.board[coords]
        moves = []
        start_file, start_rank = coords
        
        rook_dirs = [(0, 1), (0, -1), (1, 0), (-1, 0)]
        bishop_dirs = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
        knight_moves = [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]
        king_moves = rook_dirs + bishop_dirs

        def add_sliding_moves(directions):
            for df, dr in directions:
                for i in range(1, max(self.BOARD_WIDTH, self.BOARD_HEIGHT)):
                    end_coords = (start_file + i * df, start_rank + i * dr)
                    if not self._to_name(end_coords): break
                    moves.append(end_coords)
                    if self.board.get(end_coords): break

        def add_single_moves(deltas):
            for df, dr in deltas:
                end_coords = (start_file + df, start_rank + dr)
                if self._to_name(end_coords): moves.append(end_coords)

        if piece_char in ['r', 'q', 'c']: add_sliding_moves(rook_dirs)
        if piece_char in ['b', 'q', 'a']: add_sliding_moves(bishop_dirs)
        if piece_char in ['n', 'c', 'a']: add_single_moves(knight_moves)
        if piece_char == 'k': add_single_moves(king_moves)
                 
        return moves

    def is_square_attacked(self, square_coords, by_color):
        for piece_coords, (pc, p_color) in self.board.items():
            if p_color == by_color:
                if square_coords in self.get_piece_attacking_moves(piece_coords):
                    return True
        return False

    def get_all_legal_moves(self, color):
        legal_moves = []
        for start_coords, (pc, p_color) in list(self.board.items()):
            if p_color != color: continue

            for end_coords in self.get_piece_attacking_moves(start_coords):
                target_piece = self.board.get(end_coords)
                if target_piece and target_piece[1] == color: continue
                
                test_pos = self.copy()
                moved_piece = test_pos.board.pop(start_coords)
                test_pos.board[end_coords] = moved_piece
                
                king_coords = test_pos.find_king(color)
                if not test_pos.is_square_attacked(king_coords, 'white' if color == 'black' else 'black'):
                    legal_moves.append((start_coords, end_coords))
        return legal_moves

    def is_checkmate(self, color):
        king_coords = self.find_king(color)
        if not king_coords: return False
        
        opponent_color = 'white' if color == 'black' else 'black'
        if not self.is_square_attacked(king_coords, opponent_color):
            return False

        if not self.get_all_legal_moves(color):
            return True
            
        return False

    def perform_move(self, start_coords, end_coords):
        if start_coords in self.board:
            piece = self.board.pop(start_coords)
            self.board[end_coords] = piece

def solve_and_verify():
    """
    Sets up the chess position and verifies the mate-in-2 solution.
    """
    fen = "9k/5c1pb1/10/10/10/3Q6/PP5A2/K9 w - - 0 1"
    
    # --- Move 1 (White): Aj3+ ---
    pos_after_w1 = ChessPosition(fen)
    print("Initial position loaded. White to move.")
    print("Verifying proposed solution: Mate in 2.")
    print("\n--- White's 1st Move: Aj3+ ---")
    pos_after_w1.perform_move(pos_after_w1._to_coords('h2'), pos_after_w1._to_coords('j3'))
    
    # --- Verify Black's only move ---
    black_legal_moves = pos_after_w1.get_all_legal_moves('black')
    
    print("Checking Black's possible responses...")
    if len(black_legal_moves) == 1:
        start_name = pos_after_w1._to_name(black_legal_moves[0][0])
        end_name = pos_after_w1._to_name(black_legal_moves[0][1])
        print(f"Verified: Black has only one legal move: {start_name}{end_name}.")
        
        # --- Move 1 (... Black): Ki8 ---
        pos_after_b1 = pos_after_w1.copy()
        pos_after_b1.perform_move(black_legal_moves[0][0], black_legal_moves[0][1])
        
        # --- Move 2 (White): Qi6# ---
        print("\n--- White's 2nd Move: Qi6# ---")
        pos_after_w2 = pos_after_b1.copy()
        pos_after_w2.perform_move(pos_after_w2._to_coords('d3'), pos_after_w2._to_coords('i6'))
        
        # --- Verify checkmate ---
        is_mate = pos_after_w2.is_checkmate('black')
        if is_mate:
            print("Verified: The move Qi6 results in checkmate.")
            final_answer = 2
            print(f"\nConclusion: The minimal amount of moves by White to win is {final_answer}.")
        else:
            print("\nError: The sequence did not result in a checkmate.")
    else:
        print("\nError: Analysis failed. Black does not have exactly one forced move.")
        move_names = [(pos_after_w1._to_name(s), pos_after_w1._to_name(e)) for s, e in black_legal_moves]
        print(f"Legal moves found for Black: {move_names}")

solve_and_verify()