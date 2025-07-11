import pycuber as pc
import collections

def solve_f2l():
    """
    This function sets up the cube, runs a search for the shortest F2L solution,
    and prints the result.
    """

    # 1. Setup the cube and apply the scramble
    scramble_str = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    
    # pycuber requires spaces between all moves.
    scramble_str = scramble_str.replace("B'F'", "B' F'")

    cube = pc.Cube()
    cube(scramble_str)

    # 2. Reorient the cube: from (White top, Green front) to (Yellow top, Orange front)
    # This corresponds to an x2 rotation followed by a y' rotation.
    cube.perform_rotation("x2")
    cube.perform_rotation("y'")
    
    # After reorientation, the faces are:
    # U: Yellow, F: Orange, D: White, B: Red, R: Blue, L: Green

    # 3. Define functions to check if each of the 4 F2L pairs is solved
    # A pair is solved if its pieces are in the correct slot and orientation.
    
    # White-Green-Orange pair, goes in Down-Front-Left slot
    def is_wgo_solved(c):
        corner = c.get_piece("D", "F", "L")
        edge = c.get_piece("F", "L")
        return (
            corner.has_colours('W', 'G', 'O') and edge.has_colours('G', 'O') and
            corner.get_face_colour('D') == 'W' and corner.get_face_colour('F') == 'O' and corner.get_face_colour('L') == 'G' and
            edge.get_face_colour('F') == 'O' and edge.get_face_colour('L') == 'G'
        )

    # White-Green-Red pair, goes in Down-Back-Left slot
    def is_wgr_solved(c):
        corner = c.get_piece("D", "B", "L")
        edge = c.get_piece("B", "L")
        return (
            corner.has_colours('W', 'G', 'R') and edge.has_colours('G', 'R') and
            corner.get_face_colour('D') == 'W' and corner.get_face_colour('B') == 'R' and corner.get_face_colour('L') == 'G' and
            edge.get_face_colour('B') == 'R' and edge.get_face_colour('L') == 'G'
        )

    # White-Blue-Red pair, goes in Down-Back-Right slot
    def is_wbr_solved(c):
        corner = c.get_piece("D", "B", "R")
        edge = c.get_piece("B", "R")
        return (
            corner.has_colours('W', 'B', 'R') and edge.has_colours('B', 'R') and
            corner.get_face_colour('D') == 'W' and corner.get_face_colour('B') == 'R' and corner.get_face_colour('R') == 'B' and
            edge.get_face_colour('B') == 'R' and edge.get_face_colour('R') == 'B'
        )

    # White-Blue-Orange pair, goes in Down-Front-Right slot
    def is_wbo_solved(c):
        corner = c.get_piece("D", "F", "R")
        edge = c.get_piece("F", "R")
        return (
            corner.has_colours('W', 'B', 'O') and edge.has_colours('B', 'O') and
            corner.get_face_colour('D') == 'W' and corner.get_face_colour('F') == 'O' and corner.get_face_colour('R') == 'B' and
            edge.get_face_colour('F') == 'O' and edge.get_face_colour('R') == 'B'
        )

    # 4. Define goal checks for the 6 combinations of two pairs
    pair_checkers = {
        "WGO_WGR": lambda c: is_wgo_solved(c) and is_wgr_solved(c),
        "WGO_WBR": lambda c: is_wgo_solved(c) and is_wbr_solved(c),
        "WGO_WBO": lambda c: is_wgo_solved(c) and is_wbo_solved(c),
        "WGR_WBR": lambda c: is_wgr_solved(c) and is_wbr_solved(c),
        "WGR_WBO": lambda c: is_wgr_solved(c) and is_wbo_solved(c),
        "WBR_WBO": lambda c: is_wbr_solved(c) and is_wbo_solved(c),
    }

    all_moves = [pc.Move(m) for m in "U U' U2 D D' D2 L L' L2 R R' R2 F F' F2 B B' B2".split()]

    # 5. Define the BFS search function
    def bfs_solve(start_cube, goal_func):
        q = collections.deque([(start_cube.copy(), [])])
        # Using cube.formula() as the key is faster than str(cube)
        visited = {start_cube.formula()}
        
        while q:
            current_cube, path = q.popleft()
            
            if len(path) > 7: # Optimization: A solution is unlikely to be very long.
                continue

            for move in all_moves:
                if path and move.name[0] == path[-1].name[0]:
                    continue

                next_cube = current_cube.copy()(move)
                next_formula = next_cube.formula()

                if next_formula not in visited:
                    new_path = path + [move]
                    if goal_func(next_cube):
                        return new_path
                    
                    visited.add(next_formula)
                    q.append((next_cube, new_path))
        return None

    # 6. Run the search for all 6 pair combinations
    solutions = {}
    print("Searching for the shortest sequence to solve two F2L pairs...")
    for name, checker in pair_checkers.items():
        print(f"Checking combination: {name.replace('_', ' and ')}...")
        solution = bfs_solve(cube, checker)
        if solution:
            solutions[name] = solution
            print(f"Found solution of length {len(solution)}: {' '.join(str(m) for m in solution)}")

    # 7. Find and print the best solution
    if not solutions:
        print("No solution found within the search depth limit.")
        return

    best_pair_name = min(solutions, key=lambda k: len(solutions[k]))
    best_solution = solutions[best_pair_name]
    
    print("\n---")
    print(f"The best combination to solve is: {best_pair_name.replace('_', ' and ')}")
    print(f"The minimum number of moves is: {len(best_solution)}")
    print(f"The move sequence is: {' '.join(str(m) for m in best_solution)}")
    print("---")
    
    # Final answer formatting
    final_output_moves = []
    for move in best_solution:
        final_output_moves.append(str(move))
    print(f"Final equation: {' + '.join(final_output_moves)}")

    return len(best_solution)

# Execute the main function
if __name__ == '__main__':
    min_len = solve_f2l()
    print(f"<<<{min_len}>>>")
