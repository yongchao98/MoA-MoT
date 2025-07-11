import sys
import pulp

def solve_unicorn_covering():
    """
    This function formulates and solves the ILP problem to find the minimum
    number of unicorns to attack all black squares on an 8x8x8 board.
    """
    BOARD_SIZE = 8

    # Using coordinates from 1 to 8 as in standard chess notation
    all_squares = [(x, y, z) for x in range(1, BOARD_SIZE + 1)
                   for y in range(1, BOARD_SIZE + 1)
                   for z in range(1, BOARD_SIZE + 1)]

    # A square (x, y, z) is black if the sum of its coordinates is even.
    black_squares = [c for c in all_squares if sum(c) % 2 == 0]

    # Create the ILP problem
    prob = pulp.LpProblem("Unicorn_Covering_Black_Squares", pulp.LpMinimize)

    # Decision variables: u_c = 1 if a unicorn is placed on square c, 0 otherwise.
    unicorn_vars = pulp.LpVariable.dicts("Unicorn", all_squares, cat='Binary')

    # Objective function: Minimize the total number of unicorns
    prob += pulp.lpSum(unicorn_vars)

    # Constraints: Each black square must be attacked by at least one unicorn.
    # To build constraints efficiently, for each black square, we find all squares
    # from which a unicorn could attack it.
    for black_sq in black_squares:
        bx, by, bz = black_sq
        attacker_positions = []
        for s in all_squares:
            sx, sy, sz = s
            # A unicorn attacks along space diagonals. This means the absolute
            # differences in coordinates must be equal and non-zero.
            dx = abs(sx - bx)
            dy = abs(sy - by)
            dz = abs(sz - bz)
            if dx > 0 and dx == dy == dz:
                attacker_positions.append(s)
        
        # Add the constraint for this black square
        if attacker_positions:
            prob += pulp.lpSum(unicorn_vars[pos] for pos in attacker_positions) >= 1

    # Solve the problem. Suppress solver messages for cleaner output.
    solver = pulp.PULP_CBC_CMD(msg=0)
    prob.solve(solver)

    # Print the results
    if prob.status == pulp.LpStatusOptimal:
        min_unicorns = int(pulp.value(prob.objective))
        print(f"The minimum number of unicorns required is: {min_unicorns}")
        
        placements = []
        for s in all_squares:
            if unicorn_vars[s].varValue > 0.99: # Check if the variable is set to 1
                placements.append(s)

        # To satisfy the prompt "output each number in the final equation",
        # we format the result as a conceptual equation.
        print("\nThis can be represented by the equation:")
        equation_parts = [f"Unicorn_at_{c}".replace(", ", "_").replace("(", "").replace(")", "") for c in sorted(placements)]
        print(" + ".join(equation_parts) + f" = {len(placements)}")
        
        print("\nAn example of optimal unicorn placements:")
        for p in sorted(placements):
            print(f"- Unicorn at {p}")

    else:
        print(f"Could not find the optimal solution. Status: {pulp.LpStatus[prob.status]}")

if __name__ == '__main__':
    try:
        import pulp
    except ImportError:
        print("This script requires the 'pulp' library.")
        print("Please install it using: pip install pulp")
        sys.exit(1)
        
    solve_unicorn_covering()
