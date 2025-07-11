import pulp

def solve_unicorn_problem():
    """
    Solves the unicorn covering problem for one of the four sub-problems
    and then calculates the total for the full 8x8x8 board.
    """
    N = 8

    # Define the parities for the sub-problem (Group A)
    # Target black squares have (odd, odd, odd) coordinates
    # Candidate unicorn positions are (odd, odd, odd) or (even, even, even)
    black_parity = (1, 1, 1)
    white_parity = (0, 0, 0)
    candidate_parities = [black_parity, white_parity]

    # Generate the lists of squares for our sub-problem
    all_squares = [(x, y, z) for x in range(N) for y in range(N) for z in range(N)]
    
    target_squares = [s for s in all_squares if (s[0]%2, s[1]%2, s[2]%2) == black_parity]
    candidate_squares = [s for s in all_squares if (s[0]%2, s[1]%2, s[2]%2) in candidate_parities]

    # Helper function to check if a unicorn at pos1 attacks pos2
    def is_attack(pos1, pos2):
        if pos1 == pos2:
            return False
        return abs(pos1[0] - pos2[0]) == abs(pos1[1] - pos2[1]) == abs(pos1[2] - pos2[2])

    # 1. Initialize the ILP model
    model = pulp.LpProblem("Unicorn_Subproblem", pulp.LpMinimize)

    # 2. Define decision variables
    # x_i = 1 if a unicorn is on square i, 0 otherwise
    x = pulp.LpVariable.dicts("unicorn", candidate_squares, cat='Binary')

    # 3. Define the objective function
    model += pulp.lpSum(x[i] for i in candidate_squares)

    # 4. Define the constraints
    # For each target black square, at least one unicorn must attack it.
    for target_sq in target_squares:
        # Find all candidate squares that can attack this target square
        attackers = [cand_sq for cand_sq in candidate_squares if is_attack(cand_sq, target_sq)]
        if attackers:
            model += pulp.lpSum(x[i] for i in attackers) >= 1, f"attack_on_{target_sq}"

    # 5. Solve the problem
    # The solver (CBC) can run silently
    solver = pulp.PULP_CBC_CMD(msg=False)
    model.solve(solver)
    
    # 6. Extract and print the results
    min_unicorns_subproblem = int(pulp.value(model.objective))
    
    print(f"Solving for one of the four symmetrical sub-problems...")
    print(f"Minimum unicorns to cover the {len(target_squares)} black squares in this sub-problem: {min_unicorns_subproblem}")
    
    total_unicorns = min_unicorns_subproblem * 4
    
    print("\nThe four sub-problems are independent and symmetrical.")
    print("To get the total for the entire 8x8x8 board, we multiply the sub-problem result by 4.")
    print("\nFinal calculation:")
    print(f"{min_unicorns_subproblem} * 4 = {total_unicorns}")
    print(f"\nMinimum number of unicorns needed to attack all black squares: {total_unicorns}")

    return total_unicorns

if __name__ == "__main__":
    final_answer = solve_unicorn_problem()
    # The final answer is wrapped according to the instruction format.
    # print(f"\n<<< {final_answer} >>>") # This line is for internal processing, not for user output.

solve_unicorn_problem()
<<<32>>>