import pulp
import math

def solve_packing_for_box(L, W, H):
    """
    Solves the 3D packing ILP for a box of given dimensions.
    We model the box as a grid of 2x2x2 cm cells.
    """
    # Grid dimensions in terms of 2cm cells
    Gx = int(L // 2)
    Gy = int(W // 2)
    Gz = int(H // 2)

    # If the box is too small to hold even one 2x2x2 cell, no energy can be packed.
    if Gx < 1 or Gy < 1 or Gz < 1:
        return 0, 0, 0

    # Create the ILP problem
    prob = pulp.LpProblem("PackingProblem", pulp.LpMaximize)

    # --- Decision Variables ---
    # y_pqr = 1 if a small ball (1x1x1 cell block) is in cell (p,q,r)
    y_vars = pulp.LpVariable.dicts("small_ball", (range(Gx), range(Gy), range(Gz)), cat='Binary')

    # x_ijk = 1 if a large ball (2x2x2 cell block) STARTS at the corner of cell (i,j,k)
    # The start indices are limited because the 2x2x2 block needs space.
    x_vars = pulp.LpVariable.dicts("large_ball", (range(Gx - 1), range(Gy - 1), range(Gz - 1)), cat='Binary')

    # --- Objective Function ---
    # Maximize total energy: E = 1 * n1 + 20 * n2
    # where n1 = sum(y_vars) and n2 = sum(x_vars)
    prob += pulp.lpSum(y_vars) + 20 * pulp.lpSum(x_vars), "TotalEnergy"

    # --- Constraints ---
    # Each 2x2x2 cm cell in the grid can be occupied at most once, either by
    # a single small ball or as part of one large ball.
    for p in range(Gx):
        for q in range(Gy):
            for r in range(Gz):
                # Sum of all large balls that occupy cell (p,q,r)
                # A large ball started at (i,j,k) covers the 2x2x2 block of cells
                # from (i,j,k) to (i+1,j+1,k+1). So for a given cell (p,q,r),
                # we need to sum over all large ball start positions (i,j,k) such
                # that their volume contains (p,q,r).
                overlapping_x = pulp.lpSum(
                    x_vars[i][j][k]
                    for i in range(max(0, p - 1), min(p + 1, Gx - 1))
                    for j in range(max(0, q - 1), min(q + 1, Gy - 1))
                    for k in range(max(0, r - 1), min(r + 1, Gz - 1))
                )
                prob += y_vars[p][q][r] + overlapping_x <= 1, f"Cell_Constraint_{p}_{q}_{r}"

    # --- Solve the ILP ---
    # We use the default CBC solver that comes with pulp.
    # msg=False suppresses solver output.
    prob.solve(pulp.PULP_CBC_CMD(msg=False))

    # --- Extract Results ---
    if pulp.LpStatus[prob.status] == 'Optimal':
        total_energy = pulp.value(prob.objective)
        # For PuLP, need to check if varValue is not None
        num_small = sum(v.varValue for v in y_vars.values() if v.varValue is not None)
        num_large = sum(v.varValue for v in x_vars.values() if v.varValue is not None)
        return total_energy, int(round(num_small)), int(round(num_large))
    else:
        # Return a marker for failed optimization
        return -1, -1, -1

def find_optimal_container():
    """
    Searches through valid box dimensions to find the one maximizing energy.
    """
    max_energy = 0
    best_config = {}
    
    # We assume L <= W <= H to avoid redundant permutations.
    # Max possible value for the smallest dimension L: 6*L^2 <= 1050 -> L <= 13.2
    # So we iterate L up to 13.0
    l_range = [i * 0.5 for i in range(1, int(13.0 / 0.5) + 2)]

    for l_dim in l_range:
        # Max possible W for a given L: W^2 + 2*L*W - 525 <= 0
        w_limit = -l_dim + math.sqrt(l_dim**2 + 525)
        w_range = [i * 0.5 for i in range(int(l_dim / 0.5), int(w_limit / 0.5) + 2)]

        for w_dim in w_range:
            if 2 * l_dim * w_dim > 1050:
                continue

            # Max possible H for given L, W: H <= (525 - L*W)/(L+W)
            # and we require H >= W
            h_limit = (525 - l_dim * w_dim) / (l_dim + w_dim)
            h_range = [i * 0.5 for i in range(int(w_dim / 0.5), int(h_limit / 0.5) + 2)]

            for h_dim in h_range:
                surface_area = 2 * (l_dim * w_dim + l_dim * h_dim + w_dim * h_dim)
                if surface_area <= 1050:
                    energy, n1, n2 = solve_packing_for_box(l_dim, w_dim, h_dim)
                    if energy > max_energy:
                        max_energy = energy
                        best_config = {
                            "description": f"box {l_dim}x{w_dim}x{h_dim}",
                            "n1": n1,
                            "n2": n2,
                            "energy": energy,
                            "sa": surface_area
                        }

    # Print final results
    desc = best_config["description"]
    n1 = best_config["n1"]
    n2 = best_config["n2"]
    
    print(f"[{desc}]{n1};{n2}")
    
    total_energy = n1 * 1 + n2 * 20
    # As requested, output the final equation.
    print(f"Final Energy Equation: {n1} * 1 MJ + {n2} * 20 MJ = {total_energy} MJ")
    
    # The final answer in the specified format
    print(f'<<<[{desc}]{n1};{n2}>>>')


if __name__ == '__main__':
    find_optimal_container()