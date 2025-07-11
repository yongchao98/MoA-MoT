import sys
import subprocess

# Ensure the 'pulp' library is installed for Integer Linear Programming.
try:
    import pulp
except ImportError:
    print("The 'pulp' library is not installed. Attempting to install now...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pulp"])
        import pulp
        print("Successfully installed 'pulp'.")
    except Exception as e:
        print(f"Error installing pulp: {e}")
        print("Please install it manually by running: pip install pulp")
        sys.exit(1)

def solve_unicorn_problem():
    """
    Solves the 3D unicorn covering problem using Integer Linear Programming.
    """
    N = 8

    # --- Step 1: Decompose the problem ---
    # The board is split into 4 independent partitions based on unicorn moves.
    # A unicorn at (x,y,z) only moves to (x',y',z') if (x-y)%2==(x'-y')%2 and (y-z)%2==(y'-z')%2.
    # We solve the problem for one partition and multiply the result by 4.
    # We choose the partition where (x-y) and (y-z) are both even. This corresponds to
    # x, y, and z all having the same parity.
    
    print("Generating squares for one of the four partitions...")
    
    # These are all the squares available for placing unicorns in this partition
    partition_squares = []
    # These are the black squares within this partition that must be attacked
    partition_black_squares = []

    for i in range(N):
        for j in range(N):
            for k in range(N):
                # Condition for this partition: x, y, z have the same parity
                if (i % 2) == (j % 2) and (j % 2) == (k % 2):
                    square = (i, j, k)
                    partition_squares.append(square)
                    # A square is black if the sum of its coordinates is odd.
                    if (i + j + k) % 2 != 0:
                        partition_black_squares.append(square)

    print(f"Sub-problem: Cover {len(partition_black_squares)} black squares using unicorns from a set of {len(partition_squares)} squares.\n")

    # --- Step 2: Formulate the ILP for the subproblem ---
    # Define the ILP problem
    model = pulp.LpProblem("Unicorn_Subproblem", pulp.LpMinimize)

    # Define decision variables: a variable for each potential unicorn placement
    unicorn_vars = pulp.LpVariable.dicts("Unicorn", partition_squares, cat='Binary')

    # Set the objective function: minimize the number of unicorns
    model += pulp.lpSum(unicorn_vars)

    # Set the constraints: each black square must be attacked
    for b_sq in partition_black_squares:
        bx, by, bz = b_sq
        attackers = []
        for s_sq in partition_squares:
            # A unicorn at s_sq attacks b_sq if they are on a common diagonal.
            # It cannot attack its own square.
            if s_sq == b_sq:
                continue
            
            sx, sy, sz = s_sq
            if abs(sx - bx) == abs(sy - by) and abs(sy - by) == abs(sz - bz):
                attackers.append(unicorn_vars[s_sq])
        
        # Add constraint that the sum of unicorns attacking this square must be at least 1
        if attackers:
            model += pulp.lpSum(attackers) >= 1, f"Cover_Black_Square_{b_sq}"

    # --- Step 3: Solve the ILP and find the total ---
    print("Solving the subproblem with Integer Linear Programming...")
    # Use a solver (CBC is PuLP's default). msg=0 suppresses solver output.
    model.solve(pulp.PULP_CBC_CMD(msg=0))
    print("-" * 30)

    status = pulp.LpStatus[model.status]
    if status == 'Optimal':
        subproblem_solution = int(pulp.value(model.objective))
        total_unicorns = 4 * subproblem_solution

        print(f"Status for Subproblem: {status}")
        print(f"Minimum unicorns for one partition: {subproblem_solution}")

        print("\nUnicorn locations for this partition's solution:")
        solution_placements = []
        for s_sq in partition_squares:
            if unicorn_vars[s_sq].varValue > 0.5: # If unicorn is placed here
                solution_placements.append(str(s_sq))
        print(", ".join(solution_placements))
        
        # Fulfilling the request to "output each number in the final equation"
        equation = " + ".join(["1"] * subproblem_solution)
        print(f"\nEquation for the subproblem: {equation} = {subproblem_solution}")
        
        print("\n" + "-"*30)
        print(f"Since there are 4 identical partitions, the total number of unicorns is 4 * {subproblem_solution}.")
        print(f"Final Answer: The minimum number of unicorns needed is {total_unicorns}.")

        return total_unicorns
    else:
        print(f"Could not find an optimal solution. Status: {status}")
        return None

if __name__ == '__main__':
    min_unicorns = solve_unicorn_problem()
    if min_unicorns is not None:
        print(f"\n<<< {min_unicorns} >>>")
