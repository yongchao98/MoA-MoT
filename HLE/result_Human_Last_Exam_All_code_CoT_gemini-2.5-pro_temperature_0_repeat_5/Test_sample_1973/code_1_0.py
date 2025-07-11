import pulp

def solve_raumschach_unicorn_problem():
    """
    This function formulates and solves the problem of finding the minimum number of
    unicorns to attack all black squares on an 8x8x8 Raumschach board using
    Integer Linear Programming.
    """
    
    # --- Step 1: Model and solve the 4x4x4 subproblem ---
    # The overall problem decomposes into 4 independent problems on a 4x4x4 grid.
    # We solve for one 4x4x4 grid and multiply the result by 4.
    N = 4
    
    # Create the ILP problem using the PuLP library
    prob = pulp.LpProblem("Unicorn_Covering_4x4x4", pulp.LpMinimize)

    # Create a list of all cell coordinates on the 4x4x4 board
    cells = [(i, j, k) for i in range(N) for j in range(N) for k in range(N)]

    # Create binary decision variables: x_ijk = 1 if a unicorn is on cell (i,j,k), 0 otherwise.
    x = pulp.LpVariable.dicts("x", cells, cat='Binary')

    # Set the objective function: Minimize the total number of unicorns
    prob += pulp.lpSum(x[c] for c in cells), "Total_Unicorns"

    # Set the constraints: Each cell must be attacked by at least one unicorn.
    # A unicorn at (i,j,k) attacks (r,s,t) if they are on a 3D diagonal,
    # which means |i-r| == |j-s| == |k-t|. This includes the case where the
    # unicorn is on the cell itself (i=r, j=s, k=t).
    for r, s, t in cells:
        target_cell = (r, s, t)
        # Find all cells from which a unicorn can attack the target cell
        attacking_cells = [x[i, j, k] for i, j, k in cells if abs(i - r) == abs(j - s) == abs(k - t)]
        
        # Add the constraint for the target cell: the sum of unicorns on attacking
        # squares must be at least 1.
        prob += pulp.lpSum(attacking_cells) >= 1, f"Constraint_for_cell_{r}_{s}_{t}"

    # Solve the ILP problem. The solver will find the minimum value for the objective function.
    # We suppress the solver's own output for a cleaner result.
    prob.solve(pulp.PULP_CBC_CMD(msg=False))

    # Extract the result for the 4x4x4 board
    min_unicorns_4x4x4 = int(pulp.value(prob.objective))

    # --- Step 2: Calculate and print the final answer for the 8x8x8 board ---
    total_unicorns = 4 * min_unicorns_4x4x4

    print("Step-by-step explanation:")
    print("1. The problem is to find the minimum number of unicorns to attack all 256 black squares on an 8x8x8 board.")
    print("2. A key insight is that the board can be decomposed. A unicorn always stays on squares of the same color, and the set of black squares can be partitioned into 4 independent groups.")
    print("3. Each group corresponds to a unicorn covering problem on a 4x4x4 board.")
    print("4. We solve for one such 4x4x4 board using Integer Linear Programming.")
    print(f"5. The ILP solver found that the minimum number of unicorns to cover a 4x4x4 board is: {min_unicorns_4x4x4}")
    print("\nFinal Calculation:")
    print("To cover all black squares on the 8x8x8 board, we need to solve the problem for each of the 4 independent groups and sum the results.")
    print("Total Unicorns = (Unicorns for Group 1) + (Unicorns for Group 2) + (Unicorns for Group 3) + (Unicorns for Group 4)")
    print(f"The final equation is: 4 * {min_unicorns_4x4x4} = {total_unicorns}")

    print(f"\nTherefore, the minimum number of unicorns needed to attack all black squares is {total_unicorns}.")
    print(f"<<<{total_unicorns}>>>")

if __name__ == '__main__':
    # To run this code, you need to have the 'pulp' library installed.
    # You can install it using pip: pip install pulp
    solve_raumschach_unicorn_problem()