import pulp

def solve_subproblem(N, black_cell_parity, pulp_solver):
    """
    Solves the unicorn covering problem for one of the four disconnected sets of squares.

    Args:
        N (int): The size of the board (e.g., 8 for an 8x8x8 board).
        black_cell_parity (tuple): A tuple of parities, e.g., (1, 1, 1) for (odd, odd, odd).
        pulp_solver: An instance of a PuLP solver.

    Returns:
        tuple: A tuple containing the minimum number of unicorns for the subproblem
               and a list of their placement coordinates.
    """
    # p1 is the parity of the black cells to cover
    p1 = black_cell_parity
    # p2 is the parity of the white cells in the same disconnected set
    p2 = tuple((1 - p for p in p1))

    # Identify the specific black cells to be covered in this subproblem
    cells_to_cover = [
        (x, y, z) for x in range(N) for y in range(N) for z in range(N)
        if (x % 2, y % 2, z % 2) == p1
    ]

    # Identify the potential placement locations for unicorns for this subproblem
    placement_cells = [
        (x, y, z) for x in range(N) for y in range(N) for z in range(N)
        if (x % 2, y % 2, z % 2) == p1 or (x % 2, y % 2, z % 2) == p2
    ]

    # Create the ILP model
    model = pulp.LpProblem(f"Unicorn_Subproblem_{p1}", pulp.LpMinimize)

    # Define binary decision variables for each potential placement
    placement_vars = pulp.LpVariable.dicts("place", placement_cells, cat='Binary')

    # Set the objective function: minimize the number of unicorns
    model += pulp.lpSum(placement_vars[p] for p in placement_cells)

    # Add constraints: each black cell must be attacked by at least one unicorn
    for b_cell in cells_to_cover:
        bx, by, bz = b_cell
        attackers = []
        # A unicorn at p_cell attacks b_cell if they are on a common space diagonal.
        for p_cell in placement_cells:
            px, py, pz = p_cell
            if abs(px - bx) == abs(py - by) and abs(py - by) == abs(pz - bz):
                attackers.append(p_cell)
        
        # The sum of variables for cells that attack b_cell must be at least 1
        model += pulp.lpSum(placement_vars[p] for p in attackers) >= 1, f"Cover_{b_cell}"

    # Solve the ILP
    model.solve(pulp_solver)

    # Extract the results for this subproblem
    min_unicorns = int(pulp.value(model.objective))
    placements = [cell for cell in placement_cells if placement_vars[cell].value() == 1]

    return min_unicorns, placements

def solve_unicorn_problem():
    """
    Main function to solve the Raumschach unicorn problem and print the results.
    """
    # Board size
    N = 8

    # The parity patterns for the four types of black cells (where sum of components is odd)
    black_cell_parities = [(1, 1, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1)]

    total_unicorns_needed = 0
    all_placements = []
    
    # Use the default CBC solver that comes with PuLP. Suppress verbose output.
    solver = pulp.PULP_CBC_CMD(msg=False)

    # Solve each subproblem
    for parity in black_cell_parities:
        num_unicorns, placements = solve_subproblem(N, parity, solver)
        total_unicorns_needed += num_unicorns
        all_placements.extend(placements)

    # Present the final answer
    print(f"The minimum number of unicorns needed to attack all black squares is: {total_unicorns_needed}")

    print("\nOne possible set of placements for the unicorns is:")
    # Sort for consistent output
    all_placements.sort()
    for p in all_placements:
        print(f"Unicorn at {p}")
        
    print("\nThe final equation representing the sum of unicorns is:")
    # The prompt requests to output each number in the final equation.
    # This is interpreted as showing the sum of 1s for each placed unicorn.
    equation_str = " + ".join(["1"] * total_unicorns_needed)
    print(f"{equation_str} = {total_unicorns_needed}")
    
    # Finally, print the answer in the requested format
    print(f"\n<<<{total_unicorns_needed}>>>")

if __name__ == "__main__":
    solve_unicorn_problem()