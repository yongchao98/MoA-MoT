# This script uses the 'pulp' library. If you don't have it installed, run:
# pip install pulp

import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares
    on an 8x8x8 3D chessboard using Integer Linear Programming.
    """
    N = 8

    # 1. Define all squares and identify the black ones
    all_squares = [(x, y, z) for x in range(N) for y in range(N) for z in range(N)]
    black_squares = [s for s in all_squares if (s[0] + s[1] + s[2]) % 2 != 0]

    # 2. Create the ILP problem
    prob = pulp.LpProblem("Min_Unicorns_To_Cover_Black_Squares", pulp.LpMinimize)

    # 3. Define decision variables: one for each square on the board
    # The variable is 1 if a unicorn is placed, 0 otherwise.
    unicorn_vars = pulp.LpVariable.dicts("Unicorn", all_squares, cat='Binary')

    # 4. Define the objective function: minimize the total number of unicorns
    prob += pulp.lpSum([unicorn_vars[s] for s in all_squares]), "Total_Unicorns"

    # Helper function to determine if two squares are on a unicorn's path
    def can_attack(s1, s2):
        """
        Checks if a unicorn at square s1 can attack square s2.
        This is true if the squares lie on the same 3D diagonal.
        A square also "attacks" itself for covering purposes.
        """
        dx = abs(s1[0] - s2[0])
        dy = abs(s1[1] - s2[1])
        dz = abs(s1[2] - s2[2])
        return dx == dy and dx == dz

    # 5. Define the constraints: each black square must be attacked
    print("Setting up constraints...")
    for j_square in black_squares:
        # For each black square 'j', the sum of unicorns on squares 'i'
        # that can attack 'j' must be at least 1.
        covering_unicorns = [unicorn_vars[i_square] for i_square in all_squares if can_attack(i_square, j_square)]
        prob += pulp.lpSum(covering_unicorns) >= 1, f"Coverage_for_{j_square}"

    # 6. Solve the ILP problem
    print("Solving the ILP problem... (This may take a few minutes)")
    prob.solve()
    print("Solver finished.")

    # 7. Print the results
    print("\n" + "="*40)
    print("            Solution            ")
    print("="*40)
    status = pulp.LpStatus[prob.status]
    min_unicorns = int(pulp.value(prob.objective))
    
    print(f"Solver Status: {status}")
    print(f"The minimum number of unicorns needed is: {min_unicorns}")
    
    # Fulfill the request to show the final equation
    if status == 'Optimal':
        equation_numbers = ['1'] * min_unicorns
        equation_str = " + ".join(equation_numbers)
        print(f"The final equation is: {equation_str} = {min_unicorns}")

        print("\nOptimal positions for the unicorns:")
        placements = [s for s in all_squares if unicorn_vars[s].varValue == 1]
        for pos in sorted(placements):
            color = "Black" if sum(pos) % 2 != 0 else "White"
            print(f"  - {pos} (a {color} square)")
    else:
        print("Could not find the optimal solution.")
    print("="*40)


if __name__ == '__main__':
    solve_unicorn_problem()
<<<16>>>