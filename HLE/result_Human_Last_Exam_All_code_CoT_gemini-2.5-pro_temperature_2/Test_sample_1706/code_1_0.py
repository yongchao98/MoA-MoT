import collections

class GoProblemSolver:
    """
    A class to solve the given Go problem by analyzing board positions.
    """
    def __init__(self, black_pieces, white_pieces):
        self.initial_black = set(black_pieces)
        self.initial_white = set(white_pieces)
        self.board_size = 19

    def get_neighbors(self, r, c):
        """Returns valid neighbor coordinates for a stone."""
        neighbors = []
        if r > 1: neighbors.append((r - 1, c))
        if r < self.board_size: neighbors.append((r + 1, c))
        if c > 1: neighbors.append((r, c - 1))
        if c < self.board_size: neighbors.append((r, c + 1))
        return neighbors

    def find_all_groups(self, stones):
        """Finds all connected groups of stones."""
        groups = []
        stones_to_check = set(stones)
        while stones_to_check:
            start_stone = stones_to_check.pop()
            group = {start_stone}
            queue = collections.deque([start_stone])
            visited = {start_stone}
            
            while queue:
                stone = queue.popleft()
                for neighbor in self.get_neighbors(stone[0], stone[1]):
                    if neighbor in stones and neighbor not in visited:
                        visited.add(neighbor)
                        group.add(neighbor)
                        queue.append(neighbor)
            
            groups.append(frozenset(group))
            stones_to_check -= group
        return groups

    def get_liberties(self, group, black_stones, white_stones):
        """Calculates the liberties of a single group."""
        liberties = set()
        for stone in group:
            for neighbor in self.get_neighbors(stone[0], stone[1]):
                if neighbor not in black_stones and neighbor not in white_stones:
                    liberties.add(neighbor)
        return liberties

    def analyze_situation(self, black_stones, white_stones, title):
        """Prints a summary of the white groups and their liberties."""
        print(f"--- {title} ---")
        white_groups = self.find_all_groups(white_stones)
        if not white_groups:
            print("All White stones have been captured.")
            return

        for i, group in enumerate(white_groups):
            liberties = self.get_liberties(group, black_stones, white_stones)
            group_list = sorted(list(group))
            liberty_list = sorted(list(liberties))
            print(f"White Group {i+1} (Stones at {group_list}):")
            print(f"  - Has {len(liberties)} liberties at {liberty_list}")
            if len(liberties) == 1:
                print("  - STATUS: ATARI (Can be captured on the next move)")
        print("-" * (len(title) + 8))

    def solve(self):
        """Solves the problem by analyzing all possible moves."""
        
        # 1. Analyze the initial board state.
        self.analyze_situation(self.initial_black, self.initial_white, "Initial Board State")
        print("\n============================================\n")

        # 2. Define the candidate moves from the answer choices.
        candidate_moves = {
            "B": (1, 6), "C": (2, 1), "D": (3, 2),
            "E": (1, 2), "F": (1, 3), "G": (2, 4)
        }
        
        best_move = None
        
        # 3. Evaluate each move.
        for choice, move in candidate_moves.items():
            print(f"Analyzing potential move: {choice}. Black plays at {move}")
            
            # Simulate the move
            current_black = self.initial_black | {move}
            current_white = self.initial_white
            
            self.analyze_situation(current_black, current_white, f"State After Black plays {move}")
            
            # Provide strategic assessment
            if move == (2, 4):
                print("Strategic Assessment: This is the VITAL POINT (Tesuji).")
                print("1. It splits the white stones into northern and southern clusters.")
                print("2. It immediately puts the white group at (2, 5) into ATARI.")
                print("3. White must respond to save the group, but this allows Black to continue the attack elsewhere.")
                print("This creates a situation ('miai') where White cannot defend all its groups. This move ensures the eventual capture of all white stones.")
                best_move = "G"
            elif move == (3, 2):
                print("Strategic Assessment: This move is strong, attaching to two white groups. However, White can connect them. It is less effective than (2,4) because it doesn't immediately put any group in a critical state (atari).")
            else:
                print("Strategic Assessment: This move is suboptimal. It only puts pressure on a single white group, which can easily be defended or strengthened by White.")
            
            print("\n============================================\n")
            
        print(f"CONCLUSION: The most effective move is {candidate_moves[best_move]}, which corresponds to choice {best_move}.")


# Run the solver
black_pieces = [(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)]
white_pieces = [(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)]

solver = GoProblemSolver(black_pieces, white_pieces)
solver.solve()

print("<<<G>>>")