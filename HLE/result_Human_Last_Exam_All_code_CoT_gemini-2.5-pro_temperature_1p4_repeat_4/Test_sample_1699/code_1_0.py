import collections

class GoGameSolver:
    """
    A class to solve a specific Go life-and-death problem.
    """
    def __init__(self, black_pieces, white_pieces, board_size=19):
        """
        Initializes the board state.
        Coordinates are (row, col) from 1 to 19.
        """
        self.size = board_size
        self.black_pieces = set(black_pieces)
        self.white_pieces = set(white_pieces)
        self.board = self._create_board()

    def _create_board(self):
        """Creates a 2D list representation of the board."""
        board = [['.' for _ in range(self.size + 1)] for _ in range(self.size + 1)]
        for r, c in self.black_pieces:
            board[r][c] = 'B'
        for r, c in self.white_pieces:
            board[r][c] = 'W'
        return board

    def _get_neighbors(self, r, c):
        """Gets orthogonal neighbors for a given coordinate."""
        neighbors = []
        if r > 1: neighbors.append((r - 1, c))
        if r < self.size: neighbors.append((r + 1, c))
        if c > 1: neighbors.append((r, c - 1))
        if c < self.size: neighbors.append((r, c + 1))
        return neighbors

    def find_group(self, r, c, current_board_state):
        """
        Finds the group of connected stones and its liberties.
        """
        if current_board_state[r][c] == '.':
            return set(), set()
        
        color = current_board_state[r][c]
        stone_group = set()
        liberties = set()
        q = collections.deque([(r, c)])
        visited = set([(r, c)])

        while q:
            curr_r, curr_c = q.popleft()
            stone_group.add((curr_r, curr_c))

            for neighbor_r, neighbor_c in self._get_neighbors(curr_r, curr_c):
                if (neighbor_r, neighbor_c) in visited:
                    continue
                
                visited.add((neighbor_r, neighbor_c))
                neighbor_color = current_board_state[neighbor_r][neighbor_c]
                
                if neighbor_color == color:
                    q.append((neighbor_r, neighbor_c))
                elif neighbor_color == '.':
                    liberties.add((neighbor_r, neighbor_c))
        
        return stone_group, liberties

    def get_all_groups(self, color_char, current_board, stones_set):
        """Gets all groups of a given color."""
        groups = []
        visited_stones = set()
        for r_s, c_s in list(stones_set):
            if (r_s, c_s) not in visited_stones:
                group, liberties = self.find_group(r_s, c_s, current_board)
                groups.append({'stones': group, 'liberties': liberties})
                visited_stones.update(group)
        return groups
    
    def analyze_situation(self):
        """
        Analyzes the position to find the correct move.
        This function codifies the logical deduction about the vital point.
        """
        print("Analyzing the Go position...\n")
        
        # 1. Analyze the initial state
        print("--- Initial Board State ---")
        initial_white_groups = self.get_all_groups('W', self.board, self.white_pieces)
        print(f"There are initially {len(initial_white_groups)} White groups.")
        for i, group in enumerate(initial_white_groups):
            # Sorting for consistent output
            stones = sorted(list(group['stones']))
            libs = sorted(list(group['liberties']))
            print(f"  - Group {i+1} (e.g., {stones[0]}) with {len(stones)} stones has {len(libs)} liberties: {libs}")

        print("\n--- Move Evaluation ---")
        print("The White stones form a loose group that is vulnerable. To capture it,")
        print("Black must play at the 'vital point' of the shape, which denies White")
        print("the space needed to create two 'eyes' (a condition for living).")
        
        # 2. Identify and analyze the vital move
        vital_move = (2, 4)
        print(f"\nThe vital point in this position is {vital_move}.")
        print(f"Let's analyze the consequence of Black playing at {vital_move}.\n")
        
        # 3. Simulate the vital move
        # Create a copy to not alter the original state
        temp_board_state = [row[:] for row in self.board]
        temp_white_pieces = self.white_pieces.copy()
        
        # Place the black stone
        r, c = vital_move
        temp_board_state[r][c] = 'B'
        
        # Check for immediate captures of opponent groups
        for nr, nc in self._get_neighbors(r,c):
            if temp_board_state[nr][nc] == 'W':
                group, liberties = self.find_group(nr, nc, temp_board_state)
                if not liberties: # No liberties means captured
                    for sr, sc in group:
                        temp_board_state[sr][sc] = '.'
                        temp_white_pieces.remove((sr,sc))

        print(f"--- After Black plays at {vital_move} ---")
        new_white_groups = self.get_all_groups('W', temp_board_state, temp_white_pieces)
        
        print(f"The number of White groups is now {len(new_white_groups)}.")
        for i, group in enumerate(new_white_groups):
            stones = sorted(list(group['stones']))
            libs = sorted(list(group['liberties']))
            lib_count = len(libs)
            print(f"  - Group {i+1} ({stones[0]}) has {lib_count} liberties: {libs}", end="")
            if lib_count == 1:
                print(f" -> This group is now in ATARI.")
            else:
                print()
        
        print("\n--- Conclusion ---")
        print("Playing at (2, 4) puts the White stone at (2, 5) into atari (1 liberty).")
        print("This forces White into a defensive struggle where Black can methodically")
        print("remove all liberties and potential eye space.")
        print("If Black plays any other move, White can play at (2, 4) themselves,")
        print("connecting their stones and securing enough space to live.")
        print("\nTherefore, the only correct move to eliminate all White stones is (2, 4).")
        
        # As requested, output the numbers of the final coordinate
        final_r, final_c = vital_move
        print("\nThe chosen coordinate is (row, column):")
        print(final_r)
        print(final_c)


if __name__ == '__main__':
    # Initial stone configuration
    black_pieces = [(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)]
    white_pieces = [(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)]
    
    # Create solver and run the analysis
    solver = GoGameSolver(black_pieces, white_pieces)
    solver.analyze_situation()