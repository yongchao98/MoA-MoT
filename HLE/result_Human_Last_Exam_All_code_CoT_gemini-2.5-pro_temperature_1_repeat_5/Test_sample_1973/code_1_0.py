import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares on an 8x8x8 board
    using Integer Linear Programming.
    """
    # 1. Define the board and identify squares
    N = 8
    all_squares = [(i, j, k) for i in range(N) for j in range(N) for k in range(N)]
    black_squares = [s for s in all_squares if (s[0] + s[1] + s[2]) % 2 != 0]

    # 2. Create the ILP problem
    print("Setting up the Integer Linear Programming model...")
    prob = pulp.LpProblem("Unicorn_Cover", pulp.LpMinimize)

    # 3. Define decision variables
    # x[s] = 1 if a unicorn is on square s, 0 otherwise
    x = pulp.LpVariable.dicts("unicorn_at_square", all_squares, cat='Binary')

    # 4. Set the objective function: minimize the total number of unicorns
    prob += pulp.lpSum(x[s] for s in all_squares), "Total_Unicorns"

    # 5. Set the constraints: each black square must be attacked
    print(f"Generating constraints for {len(black_squares)} black squares...")
    for target_square in black_squares:
        it, jt, kt = target_square
        # A square is attacked if a unicorn is on it or can move to it.
        # Find all squares from which a unicorn can attack the target_square.
        attacking_squares = []
        for attacker_square in all_squares:
            ia, ja, ka = attacker_square
            if abs(ia - it) == abs(ja - jt) == abs(ka - kt):
                attacking_squares.append(attacker_square)
        
        # The sum of unicorns on squares that can attack the target must be >= 1
        prob += pulp.lpSum(x[s] for s in attacking_squares) >= 1, f"attack_constraint_for_{target_square}"

    # 6. Solve the problem
    print("Solving the problem... (this may take a minute or two)")
    prob.solve()

    # 7. Print the results
    print("\n--- Solution ---")
    status = pulp.LpStatus[prob.status]
    print(f"Status: {status}")

    if status == 'Optimal':
        min_unicorns = int(pulp.value(prob.objective))
        print(f"The minimum number of unicorns needed is: {min_unicorns}")

        # The user requested to output each number in the final equation
        equation_parts = ["1"] * min_unicorns
        equation_str = " + ".join(equation_parts)
        
        print("\nThe final equation is:")
        print(f"{equation_str} = {min_unicorns}")
        
        # Optionally, print the locations of the unicorns
        print("\nA possible placement for the unicorns is at coordinates:")
        unicorn_positions = []
        for s in all_squares:
            if x[s].value() is not None and x[s].value() > 0.5:
                unicorn_positions.append(s)
        for pos in sorted(unicorn_positions):
            print(pos)
            
        return min_unicorns
    else:
        print("Could not find the optimal solution.")
        return None

if __name__ == '__main__':
    # Ensure pulp is installed: pip install pulp
    try:
        result = solve_unicorn_problem()
        if result is not None:
            # The final answer format for the platform
            print(f"\n<<< {result} >>>")
    except ImportError:
        print("The 'pulp' library is required. Please install it by running: pip install pulp")
    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have a solver installed that pulp can use, like CBC (default), GLPK, or Gurobi.")
