def print_solution(grid):
    for row in grid:
        print(','.join(row))

class SudokuGrid:
    def __init__(self, initial_grid):
        self.grid = [row[:] for row in initial_grid]
        self.size = 7
        self.letters = set('abcdefg')
        # For each cell, maintain possible values
        self.possibilities = [[set(self.letters) if cell == '' else {cell} 
                             for cell in row] for row in self.grid]
        
    def update_possibilities(self):
        changed = True
        while changed:
            changed = False
            # Update based on rows and columns
            for i in range(self.size):
                for j in range(self.size):
                    if len(self.possibilities[i][j]) == 1:
                        val = next(iter(self.possibilities[i][j]))
                        # Remove from row
                        for k in range(self.size):
                            if k != j and val in self.possibilities[i][k]:
                                self.possibilities[i][k].discard(val)
                                changed = True
                        # Remove from column
                        for k in range(self.size):
                            if k != i and val in self.possibilities[k][j]:
                                self.possibilities[k][j].discard(val)
                                changed = True
            
            # Update diagonal constraint
            diag_values = set()
            for i in range(self.size):
                j = self.size - 1 - i
                if len(self.possibilities[i][j]) == 1:
                    diag_values.add(next(iter(self.possibilities[i][j])))
            
            if len(diag_values) == 1:
                diag_letter = next(iter(diag_values))
                for i in range(self.size):
                    j = self.size - 1 - i
                    if len(self.possibilities[i][j]) > 1:
                        self.possibilities[i][j] = {diag_letter}
                        changed = True
    
    def is_valid(self):
        # Check if any cell has no possibilities
        for i in range(self.size):
            for j in range(self.size):
                if len(self.possibilities[i][j]) == 0:
                    return False
        return True
    
    def find_best_cell(self):
        min_possibilities = float('inf')
        best_cell = None
        
        # First check diagonal cells
        for i in range(self.size):
            j = self.size - 1 - i
            if len(self.possibilities[i][j]) > 1:
                if len(self.possibilities[i][j]) < min_possibilities:
                    min_possibilities = len(self.possibilities[i][j])
                    best_cell = (i, j)
        
        # Then check other cells
        for i in range(self.size):
            for j in range(self.size):
                if i + j != self.size - 1 and len(self.possibilities[i][j]) > 1:
                    if len(self.possibilities[i][j]) < min_possibilities:
                        min_possibilities = len(self.possibilities[i][j])
                        best_cell = (i, j)
        
        return best_cell
    
    def solve(self):
        self.update_possibilities()
        if not self.is_valid():
            return False
            
        cell = self.find_best_cell()
        if not cell:
            return all(len(self.possibilities[i][j]) == 1 
                      for i in range(self.size) for j in range(self.size))
        
        i, j = cell
        original_possibilities = {(r,c): self.possibilities[r][c].copy() 
                                for r in range(self.size) for c in range(self.size)}
        
        for val in sorted(self.possibilities[i][j]):
            self.possibilities[i][j] = {val}
            if self.solve():
                return True
            # Restore possibilities
            for r in range(self.size):
                for c in range(self.size):
                    self.possibilities[r][c] = original_possibilities[(r,c)].copy()
        
        return False
    
    def get_solution(self):
        return [[next(iter(self.possibilities[i][j])) 
                for j in range(self.size)] for i in range(self.size)]

# Initialize grid
initial_grid = [
    ['', '', 'd', 'g', '', 'b', 'f'],
    ['', 'd', '', '', '', 'f', ''],
    ['d', '', 'g', '', '', 'f', ''],
    ['g', '', 'b', 'f', 'e', 'c', 'd'],
    ['a', 'b', 'f', '', 'c', 'd', ''],
    ['', '', '', 'c', 'd', 'g', ''],
    ['f', '', 'c', 'd', 'g', 'a', 'b']
]

solver = SudokuGrid(initial_grid)
if solver.solve():
    solution = solver.get_solution()
    print_solution(solution)
else:
    print("No solution found")