import operator
from itertools import permutations

class KenkenSolver:
    def __init__(self, size, cages):
        self.size = size
        self.cages = cages
        self.grid = [[0] * size for _ in range(size)]
        # Map each cell to its cage index for quick lookups
        self.cell_to_cage = {}
        for i, cage in enumerate(cages):
            for r, c in cage['cells']:
                self.cell_to_cage[(r, c)] = i

    def check_cage_constraint(self, cage_index):
        """Checks if a filled cage meets its constraint."""
        cage = self.cages[cage_index]
        target = cage['val']
        op_char = cage['op']
        
        values = [self.grid[r][c] for r, c in cage['cells']]
        
        # This check is for full cages only.
        if 0 in values:
            return True

        if op_char == '=':
            return values[0] == target
        
        if op_char == '+':
            return sum(values) == target
        
        if op_char == '*':
            prod = 1
            for v in values:
                prod *= v
            return prod == target
        
        # For '-' and '/', which are non-commutative, we check all permutations.
        # This makes the solver more general, though not strictly needed for this puzzle.
        if op_char == '-':
            for p in permutations(values):
                res = p[0]
                for i in range(1, len(p)):
                    res -= p[i]
                if res == target:
                    return True
            return False

        if op_char == '/':
            for p in permutations(values):
                res = p[0]
                for i in range(1, len(p)):
                    if p[i] == 0 or res % p[i] != 0:
                        res = -1 # Invalid division
                        break
                    res //= p[i]
                if res == target:
                    return True
            return False
            
        return False

    def find_empty_cell(self):
        """Finds the first empty cell (0)."""
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    return (r, c)
        return None

    def solve(self):
        """Solves the Kenken puzzle using backtracking."""
        empty_cell = self.find_empty_cell()
        if not empty_cell:
            return True  # Grid is full, puzzle solved
        
        r, c = empty_cell

        for num in range(1, self.size + 1):
            # Check row and column constraint
            is_valid_placement = (num not in self.grid[r] and
                                  num not in [self.grid[i][c] for i in range(self.size)])
            
            if is_valid_placement:
                self.grid[r][c] = num
                
                # Check cage constraint if the cage is now full
                cage_index = self.cell_to_cage[(r,c)]
                cage_cells = self.cages[cage_index]['cells']
                cage_values = [self.grid[r_cell][c_cell] for r_cell, c_cell in cage_cells]
                
                if 0 not in cage_values: # Cage is full
                    if not self.check_cage_constraint(cage_index):
                        self.grid[r][c] = 0 # Backtrack
                        continue

                if self.solve():
                    return True
                
                # Backtrack if the recursive call did not lead to a solution
                self.grid[r][c] = 0
        
        return False

# Define the cages for the corrected (solvable) puzzle
corrected_cages = [
    {'op': '+', 'val': 4, 'cells': [(0,0), (1,0)]},
    {'op': '*', 'val': 8, 'cells': [(0,1), (0,2)]},
    {'op': '*', 'val': 3, 'cells': [(0,3), (1,3)]},
    {'op': '*', 'val': 12, 'cells': [(1,1), (2,1)]},
    {'op': '+', 'val': 3, 'cells': [(1,2), (2,2)]},
    {'op': '*', 'val': 8, 'cells': [(2,0), (3,0), (3,1)]},
    {'op': '+', 'val': 6, 'cells': [(2,3), (3,3)]},
    {'op': '=', 'val': 3, 'cells': [(3,2)]}
]

# Initialize and solve the puzzle
solver = KenkenSolver(4, corrected_cages)
if solver.solve():
    top_row = solver.grid[0]
    # Print the final result without any extra formatting
    print("".join(map(str, top_row)))
