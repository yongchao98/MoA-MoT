import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares on an 8x8x8 board.
    """
    SIZE = 8

    # 1. Generate coordinates for all squares and identify black squares
    all_squares = []
    for x in range(1, SIZE + 1):
        for y in range(1, SIZE + 1):
            for z in range(1, SIZE + 1):
                all_squares.append((x, y, z))

    black_squares = [s for s in all_squares if (s[0] + s[1] + s[2]) % 2 == 0]

    # 2. Define a function to check for unicorn attacks
    def can_attack(s1, s2):
        """Checks if a unicorn at square s1 can attack square s2."""
        # A piece on a square "attacks" that same square for the purpose of coverage.
        if s1 == s2:
            return True
        # A unicorn moves along space diagonals, meaning the change in each coordinate is the same.
        dx = abs(s1[0] - s2[0])
        dy = abs(s1[1] - s2[1])
        dz = abs(s1[2] - s2[2])
        return dx > 0 and dx == dy and dx == dz

    # 3. Set up the Integer Linear Programming problem
    prob = pulp.LpProblem("Unicorn_Cover", pulp.LpMinimize)

    # 4. Define decision variables
    # x_s = 1 if a unicorn is on square s, 0 otherwise
    x = pulp.LpVariable.dicts("unicorn_on_square", all_squares, cat='Binary')

    # 5. Define the objective function
    # Minimize the total number of placed unicorns
    prob += pulp.lpSum(x[s] for s in all_squares)

    # 6. Define the constraints
    # For each black square, at least one unicorn must attack it.
    for b in black_squares:
        attacking_squares = [s for s in all_squares if can_attack(s, b)]
        prob += pulp.lpSum(x[s] for s in attacking_squares) >= 1, f"cover_black_square_{b}"

    # 7. Solve the ILP
    # PuLP's default CBC solver is used. It might take a moment to run.
    print("Solving... (this may take a minute)")
    status = prob.solve()

    # 8. Print the results
    print(f"\nSolver status: {pulp.LpStatus[status]}")

    if pulp.LpStatus[status] == 'Optimal':
        min_unicorns = int(pulp.value(prob.objective))
        print(f"Minimum number of unicorns needed to attack all {len(black_squares)} black squares: {min_unicorns}")

        unicorn_locations = []
        for s in all_squares:
            if x[s].varValue > 0.5: # Check if the variable is 1
                unicorn_locations.append(s)

        print("\nThe final equation is the sum of one for each unicorn placed:")
        equation_parts = ["1" for _ in unicorn_locations]
        equation = " + ".join(equation_parts) + f" = {min_unicorns}"
        print(equation)
        
        # Uncomment the following lines to see the exact locations of the unicorns
        # print("\nOptimal unicorn locations (x, y, z):")
        # for loc in sorted(unicorn_locations):
        #     print(loc)
            
        return min_unicorns
    else:
        print("An optimal solution could not be found.")
        return None

if __name__ == '__main__':
    result = solve_unicorn_problem()
    if result is not None:
        print(f"\n<<< {result} >>>")
