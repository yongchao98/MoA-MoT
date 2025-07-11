import math

def solve_grid_path_problem():
    """
    Calculates the number of unique paths from (0,0) to (4,8) with the given constraints
    using the Principle of Inclusion-Exclusion.
    """
    r_steps = 4
    u_steps = 8
    total_steps = r_steps + u_steps
    max_consecutive = 3

    # Step 1: Calculate the total number of paths without any restrictions.
    # This is C(total_steps, r_steps) or C(12, 4).
    total_paths = math.comb(total_steps, r_steps)

    print(f"The goal is to get from (0,0) to ({r_steps},{u_steps}) with no more than {max_consecutive} consecutive moves in the same direction.")
    print("-" * 20)
    print("Step 1: Calculate the total number of paths without restrictions.")
    print(f"Any path has {r_steps} right moves and {u_steps} up moves, for a total of {total_steps} steps.")
    print(f"The total number of paths is C({total_steps}, {r_steps}) = {total_paths}.")
    print("-" * 20)

    # Step 2: Calculate N(A), the number of paths with a forbidden R-run (RRRR).
    # Since there are only 4 R-moves, this means the path contains the block 'RRRR'.
    # We treat 'RRRR' as a single unit and arrange it with the 8 'U' moves.
    # This is equivalent to arranging 1 'RRRR' block and 8 'U's, which is C(1+8, 1).
    paths_with_rrrr = math.comb(1 + u_steps, 1)
    
    print("Step 2: Calculate N(A), the number of paths with a forbidden 'Right' run (RRRR).")
    print("We treat 'RRRR' as a single block and arrange it with the 8 'U' moves.")
    print(f"This means arranging 9 items (1 'RRRR' block, 8 'U's), so N(A) = C(9, 1) = {paths_with_rrrr}.")
    print("-" * 20)

    # Step 3: Calculate N(B), the number of paths with a forbidden U-run (UUUU or more).
    # We can model this by placing the 8 'U's into the 5 slots created by the 4 'R's: _ R _ R _ R _ R _
    # Let u1, u2, u3, u4, u5 be the number of U's in each slot. u1+u2+u3+u4+u5 = 8.
    # We need to find the number of solutions where at least one u_i >= 4.
    # Using inclusion-exclusion on the variables:
    # N(>=1 var is >=4) = C(5,1)*C(8-4+5-1, 5-1) = C(5,1)*C(8,4) = 5 * 70 = 350
    # N(>=2 vars are >=4) = C(5,2)*C(8-8+5-1, 5-1) = C(5,2)*C(4,4) = 10 * 1 = 10
    # N(B) = 350 - 10 = 340
    n_b_term1 = math.comb(r_steps + 1, 1) * math.comb(u_steps - (max_consecutive + 1) + (r_steps + 1) - 1, r_steps)
    n_b_term2 = math.comb(r_steps + 1, 2) * math.comb(u_steps - 2 * (max_consecutive + 1) + (r_steps + 1) - 1, r_steps)
    paths_with_uuuu = n_b_term1 - n_b_term2

    print("Step 3: Calculate N(B), the number of paths with a forbidden 'Up' run (4 or more consecutive 'U's).")
    print("This is calculated by finding arrangements of R's and U's with at least one run 'UUUU...'.")
    print(f"Using combinatorial formulas (inclusion-exclusion on slots), we find N(B) = {paths_with_uuuu}.")
    print("-" * 20)

    # Step 4: Calculate N(A and B), the number of paths with both forbidden runs.
    # A path with 'RRRR' must place the 8 'U's around this block (_ RRRR _).
    # Let u1 + u2 = 8. To avoid a 'UUUU' run, we'd need u1<=3 and u2<=3, which is impossible.
    # Therefore, any path with 'RRRR' must also have a run of at least 4 'U's.
    # This means the set of paths with 'RRRR' is a subset of the paths with 'UUUU...'.
    # So, N(A and B) = N(A).
    paths_with_both = paths_with_rrrr

    print("Step 4: Calculate N(A and B), the number of paths with both forbidden runs.")
    print("Any path containing the 'RRRR' block must have a run of 4 or more 'U's.")
    print(f"Therefore, N(A and B) = N(A) = {paths_with_both}.")
    print("-" * 20)

    # Step 5: Apply the Principle of Inclusion-Exclusion.
    # Valid Paths = Total - (N(A) + N(B) - N(A and B))
    invalid_paths = paths_with_rrrr + paths_with_uuuu - paths_with_both
    final_answer = total_paths - invalid_paths

    print("Step 5: Apply the Principle of Inclusion-Exclusion to find the final answer.")
    print("Valid Paths = Total Paths - (N(A) + N(B) - N(A and B))")
    print(f"Valid Paths = {total_paths} - ({paths_with_rrrr} + {paths_with_uuuu} - {paths_with_both})")
    print(f"Valid Paths = {total_paths} - {invalid_paths}")
    print(f"The number of unique ways is {final_answer}.")

solve_grid_path_problem()
<<<155>>>