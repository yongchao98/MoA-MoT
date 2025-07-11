import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares on an 8x8x8 board using ILP.
    """
    N = 8
    
    # Step 1: Define the squares for one of the four independent classes.
    # We choose the class where (x+y)%2 == 0 and (y+z)%2 == 0.
    # A unicorn in this class can only attack other squares in this class.
    class_squares = []
    for x in range(N):
        for y in range(N):
            for z in range(N):
                if (x + y) % 2 == 0 and (y + z) % 2 == 0:
                    class_squares.append((x, y, z))

    # Identify the black squares within this class that need to be attacked.
    # A square (x,y,z) is black if (x+y+z) is odd.
    black_squares_to_cover = [s for s in class_squares if (s[0] + s[1] + s[2]) % 2 == 1]

    # The potential positions for our unicorns are all squares in the class.
    potential_unicorn_positions = class_squares

    # Step 2: Pre-compute which positions attack which black squares.
    # This helps in building the constraints efficiently.
    attack_map = {target: [] for target in black_squares_to_cover}
    for unicorn_pos in potential_unicorn_positions:
        for target_pos in black_squares_to_cover:
            dx = abs(unicorn_pos[0] - target_pos[0])
            dy = abs(unicorn_pos[1] - target_pos[1])
            dz = abs(unicorn_pos[2] - target_pos[2])
            # Check for a valid unicorn move
            if dx == dy and dy == dz:
                attack_map[target_pos].append(unicorn_pos)

    # Step 3: Set up the Integer Linear Programming problem.
    prob = pulp.LpProblem("Unicorn_Cover_Subproblem", pulp.LpMinimize)

    # Decision variables: a binary variable for each potential unicorn position.
    # u_pos = 1 if we place a unicorn at pos, 0 otherwise.
    unicorn_vars = pulp.LpVariable.dicts("Unicorn", potential_unicorn_positions, cat='Binary')

    # Objective function: Minimize the total number of placed unicorns.
    prob += pulp.lpSum(unicorn_vars[pos] for pos in potential_unicorn_positions)

    # Constraints: Each black square must be attacked by at least one unicorn.
    for target_pos in black_squares_to_cover:
        # Sum of unicorns that can attack this target_pos must be >= 1
        prob += pulp.lpSum(unicorn_vars[attacker_pos] for attacker_pos in attack_map[target_pos]) >= 1, f"Cover_{target_pos}"

    # Step 4: Solve the ILP problem for the subproblem.
    # The solver will be called automatically by PuLP. A message may be printed to the console.
    prob.solve()

    # The minimum number of unicorns for one subproblem
    min_unicorns_subproblem = int(pulp.value(prob.objective))

    # Step 5: The total number is 4 times the result of the subproblem.
    total_unicorns = 4 * min_unicorns_subproblem

    print(f"The board can be decomposed into 4 independent subproblems.")
    print(f"Minimum unicorns needed for one subproblem: {min_unicorns_subproblem}")
    print(f"The final equation is: 4 * {min_unicorns_subproblem} = {total_unicorns}")
    print(f"Total minimum number of unicorns needed: {total_unicorns}")

solve_unicorn_problem()