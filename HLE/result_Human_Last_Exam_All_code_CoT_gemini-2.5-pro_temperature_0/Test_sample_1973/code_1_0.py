import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares on an 8x8x8 board
    using Integer Linear Programming.
    """
    # Plan:
    # 1. The 3D board can be partitioned into 4 independent sets of squares based on unicorn movement.
    #    A unicorn's move preserves the parity of the differences between coordinates.
    # 2. The problem of covering all black squares breaks down into 4 identical, independent subproblems.
    # 3. We solve one subproblem using ILP to find the minimum unicorns for one partition.
    # 4. The final answer is 4 times the result of the subproblem.

    # --- Subproblem: Cover all black squares in the first partition ---
    # A square (x, y, z) is in the first partition if (x-y) and (x-z) are both even.

    # 1. Define board parameters and identify the first partition
    N = 8
    all_squares = [(x, y, z) for x in range(N) for y in range(N) for z in range(N)]

    # Partition C00: (x-y)%2==0 and (x-z)%2==0
    partition_squares = [s for s in all_squares if (s[0] - s[1]) % 2 == 0 and (s[0] - s[2]) % 2 == 0]

    # 2. Identify black squares within this partition
    # A square is black if (x+y+z) is odd.
    black_squares_in_partition = [s for s in partition_squares if (s[0] + s[1] + s[2]) % 2 == 1]

    # 3. Pre-compute which squares attack which black squares within the partition
    attacks = {}
    for s1 in partition_squares:
        attacked_black_squares = []
        for s2 in black_squares_in_partition:
            # A unicorn at s1 attacks s2 if they are on the same diagonal,
            # or if s1 is the same square as s2 (and s1 is black).
            # A unicorn at a position "attacks" that position by occupation.
            if s1 == s2:
                attacked_black_squares.append(s2)
                continue
            
            dx = abs(s1[0] - s2[0])
            dy = abs(s1[1] - s2[1])
            dz = abs(s1[2] - s2[2])
            
            # Unicorn move condition
            if dx > 0 and dx == dy and dx == dz:
                attacked_black_squares.append(s2)
        attacks[s1] = attacked_black_squares

    # 4. Set up the ILP model for the subproblem
    prob = pulp.LpProblem("Unicorn_Subproblem", pulp.LpMinimize)

    # Decision variables: x_s = 1 if a unicorn is on square s in the partition
    var_keys = [f"{s[0]}_{s[1]}_{s[2]}" for s in partition_squares]
    x = pulp.LpVariable.dicts("x", var_keys, cat='Binary')

    # Objective function: Minimize the number of unicorns in this partition
    prob += pulp.lpSum([x[f"{s[0]}_{s[1]}_{s[2]}"] for s in partition_squares]), "Unicorns_In_Partition"

    # Constraints: Each black square in the partition must be covered by at least one unicorn
    for b in black_squares_in_partition:
        b_key = f"{b[0]}_{b[1]}_{b[2]}"
        
        # Find all squares s in the partition that can attack b
        potential_attacker_positions = [s for s in partition_squares if b in attacks[s]]
                
        constraint_expr = pulp.lpSum([x[f"{s[0]}_{s[1]}_{s[2]}"] for s in potential_attacker_positions])
        prob += constraint_expr >= 1, f"Cover_Black_Square_{b_key}"

    # 5. Solve the ILP for the subproblem
    # The solver will run quietly without printing its own logs.
    prob.solve(pulp.PULP_CBC_CMD(msg=False))

    # 6. Calculate and print the final result
    if prob.status == pulp.LpStatusOptimal:
        unicorns_for_partition = int(pulp.value(prob.objective))
        total_unicorns = 4 * unicorns_for_partition

        # The prompt asks to output each number in the final equation.
        print(f"The 8x8x8 board can be divided into 4 independent sub-problems.")
        print(f"The minimum number of unicorns to cover the black squares in one sub-problem is: {unicorns_for_partition}")
        print(f"Therefore, the total minimum number of unicorns is 4 * {unicorns_for_partition} = {total_unicorns}")
        print(f"<<<{total_unicorns}>>>")
    else:
        print("Could not find the optimal solution.")
        print(f"Solver status: {pulp.LpStatus[prob.status]}")

if __name__ == '__main__':
    solve_unicorn_problem()