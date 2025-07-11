import collections

# Global cache for area decompositions to avoid recomputing
area_decomposition_cache = {}

def find_area_decompositions(target_area, tiles):
    """
    Finds all multisets of tiles that sum up to the target_area.
    This is a variation of the change-making problem.
    Uses memoization to speed up computation.
    """
    if target_area in area_decomposition_cache:
        return area_decomposition_cache[target_area]

    if target_area == 0:
        return [collections.Counter()]
    if target_area < 0 or not tiles:
        return []

    # Unpack the first tile from the list
    tile_side, tile_area = tiles[0]
    remaining_tiles = tiles[1:]
    
    decompositions = []
    
    # Iterate on the number of times we can use the current tile
    for num_tiles in range(target_area // tile_area + 1):
        current_area = num_tiles * tile_area
        current_multiset = collections.Counter({tile_side: num_tiles})
        
        remaining_area = target_area - current_area
        
        # Recursively find decompositions for the rest of the area
        sub_decompositions = find_area_decompositions(remaining_area, remaining_tiles)
        
        for sub_deco in sub_decompositions:
            decompositions.append(current_multiset + sub_deco)
            
    area_decomposition_cache[target_area] = decompositions
    return decompositions

class TilingSolver:
    """
    Solves if a rectangle can be tiled with a given multiset of squares.
    """
    def __init__(self, L, W, tile_counts):
        self.L = L
        self.W = W
        self.grid = [[0] * W for _ in range(L)]
        self.tile_counts = tile_counts
        # Sort tiles by size, largest first, for the backtracking algorithm
        self.tile_sides = sorted(tile_counts.keys(), reverse=True)

    def solve(self):
        """Initial call to the recursive solver."""
        return self._solve_recursive(self.tile_counts)

    def find_empty_cell(self):
        """Find the top-most, left-most empty cell."""
        for r in range(self.L):
            for c in range(self.W):
                if self.grid[r][c] == 0:
                    return r, c
        return None, None
        
    def _solve_recursive(self, current_tile_counts):
        """Backtracking algorithm to fit tiles into the grid."""
        r, c = self.find_empty_cell()
        if r is None:
            # No empty cells, so the grid is fully tiled
            return True

        for s in self.tile_sides:
            if current_tile_counts[s] > 0:
                # Check if the square s x s can be placed at (r, c)
                if self.can_place(r, c, s):
                    # Place the tile
                    self.place(r, c, s, s)
                    current_tile_counts[s] -= 1
                    
                    # Recurse
                    if self._solve_recursive(current_tile_counts):
                        return True
                    
                    # Backtrack
                    self.place(r, c, s, 0)
                    current_tile_counts[s] += 1
        return False

    def can_place(self, r, c, s):
        """Check if an s x s square fits at (r, c)."""
        if r + s > self.L or c + s > self.W:
            return False
        for i in range(r, r + s):
            for j in range(c, c + s):
                if self.grid[i][j] != 0:
                    return False
        return True

    def place(self, r, c, s, fill_value):
        """Place or remove a square on the grid."""
        for i in range(r, r + s):
            for j in range(c, c + s):
                self.grid[i][j] = fill_value

def find_smallest_fault_free_tileable_rectangle():
    """Main function to find the solution."""
    
    # Tile properties: (side_length, area)
    tiles = [(7, 49), (5, 25), (3, 9), (2, 4)]
    
    # Generate candidate dimensions L, W that are not divisible by 2, 3, 5, 7
    primes = {2, 3, 5, 7}
    allowed_dims = []
    for n in range(2, 50): # Search up to a reasonable dimension limit
        if all(n % p != 0 for p in primes):
            allowed_dims.append(n)
    
    # Create candidate (L, W) pairs and sort them by area
    candidates = []
    for l in allowed_dims:
        for w in allowed_dims:
            if l <= w:
                candidates.append((l, w))
    
    candidates.sort(key=lambda x: x[0] * x[1])

    for L, W in candidates:
        area = L * W
        print(f"Checking rectangle: {L}x{W} (Area: {area})")

        # Clear cache for the new area
        global area_decomposition_cache
        area_decomposition_cache = {}
        decompositions = find_area_decompositions(area, tiles)
        
        if not decompositions:
            print(f"  No area decomposition found for {area}.")
            continue
        
        print(f"  Found {len(decompositions)} area decompositions. Checking for tilings...")
        
        for deco in decompositions:
            # We only need one successful tiling
            solver = TilingSolver(L, W, deco)
            if solver.solve():
                print(f"\nSolution found!")
                print(f"The smallest rectangle is {L}x{W}.")
                print("The area of this rectangle is the product of its sides.")
                print(f"{L} * {W} = {area}")
                return area
    
    print("No solution found within the searched candidates.")
    return None

if __name__ == "__main__":
    area = find_smallest_fault_free_tileable_rectangle()
    if area:
        print(f"\nFinal Answer: The area is {area}.")
        # The required final output format for the platform
        print(f"\n<<<{area}>>>")
