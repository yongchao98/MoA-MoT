import collections
import pycuber as pc

def count_solved_f2l_pairs(cube):
    """
    Counts the number of solved F2L pairs for the Y-top, O-front orientation.
    A pair is solved if its edge and corner are in the correct place and orientation
    relative to the white (D) face.
    """
    solved_count = 0
    # Expected center colors for this orientation:
    # U=Yellow, D=White, F=Orange, B=Red, R=Blue, L=Green

    # Pair 1: Front-Right slot (Orange-Blue edge, White-Orange-Blue corner)
    is_p1_corner_ok = (cube.D[2][2].colour == 'white' and
                       cube.F[2][2].colour == 'orange' and
                       cube.R[2][0].colour == 'blue')
    is_p1_edge_ok = (cube.F[1][2].colour == 'orange' and
                     cube.R[1][0].colour == 'blue')
    if is_p1_corner_ok and is_p1_edge_ok:
        solved_count += 1

    # Pair 2: Front-Left slot (Orange-Green edge, White-Orange-Green corner)
    is_p2_corner_ok = (cube.D[2][0].colour == 'white' and
                       cube.F[2][0].colour == 'orange' and
                       cube.L[2][2].colour == 'green')
    is_p2_edge_ok = (cube.F[1][0].colour == 'orange' and
                     cube.L[1][2].colour == 'green')
    if is_p2_corner_ok and is_p2_edge_ok:
        solved_count += 1

    # Pair 3: Back-Left slot (Red-Green edge, White-Red-Green corner)
    is_p3_corner_ok = (cube.D[0][0].colour == 'white' and
                       cube.B[2][2].colour == 'red' and
                       cube.L[2][0].colour == 'green')
    is_p3_edge_ok = (cube.B[1][2].colour == 'red' and
                     cube.L[1][0].colour == 'green')
    if is_p3_corner_ok and is_p3_edge_ok:
        solved_count += 1

    # Pair 4: Back-Right slot (Red-Blue edge, White-Red-Blue corner)
    is_p4_corner_ok = (cube.D[0][2].colour == 'white' and
                       cube.B[2][0].colour == 'red' and
                       cube.R[2][2].colour == 'blue')
    is_p4_edge_ok = (cube.B[1][0].colour == 'red' and
                     cube.R[1][2].colour == 'blue')
    if is_p4_corner_ok and is_p4_edge_ok:
        solved_count += 1

    return solved_count

def solve_f2l():
    """
    Main function to set up the cube and find the shortest F2L solution.
    """
    # 1. Initial setup
    my_cube = pc.Cube()
    scramble = pc.Formula("R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B'F' U' R L' D R' B R F2 L' F2 L D")
    my_cube(scramble)

    # 2. Re-orient to Yellow top, Orange front (x2 y')
    my_cube("x2 y'")

    # 3. Check initial state
    initial_solved_count = count_solved_f2l_pairs(my_cube)
    print(f"Initially, {initial_solved_count} F2L pair(s) are solved.")

    if initial_solved_count >= 2:
        print("A solution of length 0 is found as two pairs are already solved.")
        print("The minimum length is: 0")
        print("The sequence of moves is: (None)")
        return

    # 4. Breadth-First Search for the shortest path
    # A queue of (cube_state, path_taken)
    queue = collections.deque([(my_cube, [])])
    # A set to store visited states to avoid redundant work
    visited = {str(my_cube)}
    # Moves to try (all faces except D, since the white cross is there)
    moves_to_try = []
    for face in "RUFLB":
        moves_to_try.extend([face, face + "'", face + "2"])

    while queue:
        current_cube, path = queue.popleft()

        # Try every possible next move
        for move in moves_to_try:
            next_cube = current_cube.copy()
            next_cube(move)
            
            state_str = str(next_cube)
            if state_str in visited:
                continue
            
            new_path = path + [move]
            
            # Check if we reached the goal state
            if count_solved_f2l_pairs(next_cube) >= 2:
                solution_path_str = " ".join(new_path)
                print(f"Solution found in {len(new_path)} moves: {solution_path_str}")
                print(f"The minimum length is: {len(new_path)}")
                # Print the final "equation" as requested
                print("The sequence of moves is:")
                for m in new_path:
                    print(m, end=" ")
                print()
                return

            visited.add(state_str)
            queue.append((next_cube, new_path))

# Run the solver
solve_f2l()