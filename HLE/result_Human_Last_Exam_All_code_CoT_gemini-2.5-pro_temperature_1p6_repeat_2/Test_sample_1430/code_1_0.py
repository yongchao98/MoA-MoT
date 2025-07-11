import sys

def solve_neron_severi_rank():
    """
    Calculates the smallest and largest possible rank of the Neron-Severi group
    of X = C^(15), the 15th symmetric power of a genus 3 Riemann surface C.
    """
    # Step 1: Define problem parameters
    g = 3  # Genus of the curve C
    d = 15 # Degree of the symmetric power X = C^(d)

    print("The problem is to find the smallest and largest possible rank of the Neron-Severi group of X = C^(d).")
    print(f"Given: Genus of the curve C is g = {g}.")
    print(f"       The variety X is the d-th symmetric power, with d = {d}.")
    print("-" * 50)

    # Step 2 & 3: State the relevant formula and check its condition
    print("For a symmetric power C^(d) of a curve C of genus g, there is a formula for the rank of its Neron-Severi group (the Picard number rho):")
    print("rho(C^(d)) = rho(C) + rho(J(C))")
    print("This formula is valid if d >= 2*g - 1.")

    condition_val = 2 * g - 1
    print(f"Checking the condition: We need d >= 2*g - 1.")
    print(f"Here, d = {d} and 2*g - 1 = 2*{g} - 1 = {condition_val}.")

    if d < condition_val:
        print(f"Error: The condition {d} >= {condition_val} is not met. The formula cannot be applied.", file=sys.stderr)
        return

    print(f"Since {d} >= {condition_val}, the condition holds and the formula is applicable.")
    print("-" * 50)

    # Step 4: Determine the rank for the curve C
    rho_C = 1
    print(f"For any smooth projective curve C, the rank of its Neron-Severi group is constant.")
    print(f"rho(C) = {rho_C}.")
    print("-" * 50)

    # Step 5: Determine the range for the rank of the Jacobian J(C)
    print("The problem now depends on the possible range for rho(J(C)), the rank of the Neron-Severi group of the Jacobian of a genus 3 curve.")

    # Smallest rank for J(C)
    rho_min_J_C = 1
    print(f"\nThe smallest possible rank occurs for a 'very general' curve C, where the endomorphism ring of J(C) is just Z.")
    print(f"In this case, rho_min(J(C)) = {rho_min_J_C}.")

    # Largest rank for J(C)
    rho_max_J_C = 7
    print(f"\nThe largest possible rank is achieved for special curves whose Jacobians have maximal complex multiplication.")
    print(f"For a Jacobian of a genus 3 curve, the maximum known rank is 7.")
    print(f"So, rho_max(J(C)) = {rho_max_J_C}.")
    print("-" * 50)

    # Step 6: Calculate the final results
    print("We can now compute the smallest and largest possible ranks for rho(X).")

    # Smallest rank for X
    smallest_rank = rho_C + rho_min_J_C
    print("\nSmallest Rank Calculation:")
    print(f"rho_min(X) = rho(C) + rho_min(J(C))")
    # As requested, printing each number in the final equation
    print(f"rho_min(X) = {rho_C} + {rho_min_J_C} = {smallest_rank}")


    # Largest rank for X
    largest_rank = rho_C + rho_max_J_C
    print("\nLargest Rank Calculation:")
    print(f"rho_max(X) = rho(C) + rho_max(J(C))")
    # As requested, printing each number in the final equation
    print(f"rho_max(X) = {rho_C} + {rho_max_J_C} = {largest_rank}")
    print("-" * 50)

    print(f"Final Answer: The smallest possible rank is {smallest_rank} and the largest is {largest_rank}.")

solve_neron_severi_rank()