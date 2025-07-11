import collections
import pycuber as pc
import sys

def solve_f2l_pairs():
    """
    This script finds the shortest move sequence to solve two F2L pairs
    from a given scrambled state.
    """

    # This helper function checks the F2L state of the cube. It is written
    # for the specific orientation where Johnny picks up the cube:
    # Yellow on top, Orange on front.
    # The resulting centers are: U=Yellow, D=White, F=Orange, B=Red, R=Blue, L=Green.
    def count_solved_f2l_pairs(cube):
        solved_count = 0
        
        # Pair 1: White-Orange-Blue corner, Orange-Blue edge.
        # Home position: Corner at DFR, Edge at FR.
        try:
            c = cube.get_piece('W','O','B')
            e = cube.get_piece('O','B')
            if c.pos == 'DFR' and e.pos == 'FR':
                if c.facings['D'] == 'W' and e.facings['F'] == 'O':
                    solved_count += 1
        except KeyError:
            pass # Piece not found with these exact colors
            
        # Pair 2: White-Blue-Red corner, Blue-Red edge.
        # Home position: Corner at DRB, Edge at RB.
        try:
            c = cube.get_piece('W','B','R')
            e = cube.get_piece('B','R')
            if c.pos == 'DRB' and e.pos == 'RB':
                if c.facings['D'] == 'W' and c.facings['R'] == 'B':
                    solved_count += 1
        except KeyError:
            pass

        # Pair 3: White-Red-Green corner, Red-Green edge.
        # Home position: Corner at DBL, Edge at BL.
        try:
            c = cube.get_piece('W','R','G')
            e = cube.get_piece('R','G')
            if c.pos == 'DBL' and e.pos == 'BL':
                if c.facings['D'] == 'W' and c.facings['B'] == 'R':
                    solved_count += 1
        except KeyError:
            pass

        # Pair 4: White-Green-Orange corner, Green-Orange edge.
        # Home position: Corner at DLF, Edge at LF.
        try:
            c = cube.get_piece('W','G','O')
            e = cube.get_piece('G','O')
            if c.pos == 'DLF' and e.pos == 'LF':
                if c.facings['D'] == 'W' and c.facings['L'] == 'G':
                    solved_count += 1
        except KeyError:
            pass
        
        return solved_count

    # 1. Set up the initial state from the scramble.
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B'F' U' R L' D R' B R F2 L' F2 L D"
    cube = pc.Cube()
    cube(scramble)

    # 2. Reorient the cube to match how Johnny holds it.
    # "yellow on the top and orange on the front" requires an 'x2 y' rotation.
    cube("x2 y")
    
    # 3. Prepare for the Breadth-First Search (BFS).
    moves = ["U", "U'", "U2", "D", "D'", "D2", "L", "L'", "L2", "R", "R'", "R2", "F", "F'", "F2", "B", "B'", "B2"]
    
    # The queue stores tuples of (cube_object, path_as_list_of_moves).
    queue = collections.deque([(cube, [])])
    # The visited set stores string representations of cube states to avoid cycles.
    visited = {str(cube)}
    
    # Check if the initial state already solves the problem.
    initial_solved_pairs = count_solved_f2l_pairs(cube)
    if initial_solved_pairs >= 2:
        print("The cube already has two or more F2L pairs solved.")
        print("Minimum number of moves: 0")
        print("Sequence: (None)")
        return
        
    # 4. Execute the BFS.
    # Set a maximum search depth to prevent extremely long runtimes.
    max_depth = 6
    print("Searching for the shortest sequence...")
    while queue:
        current_cube, current_path = queue.popleft()
        
        if len(current_path) >= max_depth:
            print(f"Search exceeded max depth of {max_depth}. No solution found.")
            break
            
        for move in moves:
            next_cube = current_cube.copy()
            next_cube(move)
            
            state_str = str(next_cube)
            if state_str not in visited:
                new_path = current_path + [move]
                
                # Check if this new state is the goal.
                if count_solved_f2l_pairs(next_cube) >= 2:
                    solution_path = " ".join(new_path)
                    solution_length = len(new_path)
                    
                    print(f"\nSolution Found!")
                    print(f"The minimum length is: {solution_length}")
                    print(f"The sequence of moves is:")
                    # The prompt "output each number in the final equation" is interpreted
                    # as printing the standard move notation, which includes numbers (e.g., L2).
                    print(f"{solution_path}")
                    return

                visited.add(state_str)
                queue.append((next_cube, new_path))

if __name__ == "__main__":
    try:
        solve_f2l_pairs()
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        print("This might be due to a missing library. Please try running 'pip install pycuber'.", file=sys.stderr)
