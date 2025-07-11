import math

def solve_chair_puzzle():
    """
    Solves the chair puzzle by analyzing the rules and finding the optimal strategy.
    """
    num_chairs = 20

    print("--- Puzzle Analysis ---")
    print("The goal is to find the maximum number of occupied chairs in a row of 20.")
    print("Rule: A person sits in an empty chair. If they have neighbors, one neighbor leaves.")
    print("\nLet's analyze the change in the number of occupied chairs:")
    print("1. A person sits down: The count of occupied chairs increases by 1.")
    print("2. A neighbor stands up: The count of occupied chairs decreases by 1.")
    print("\nThere are two key scenarios for a person sitting down:")
    print("  - Scenario A: The chosen empty chair has NO occupied neighbors.")
    print("    In this case, no one needs to leave. The net change is +1.")
    print("  - Scenario B: The chosen empty chair has at least ONE occupied neighbor.")
    print("    In this case, one of the neighbors must leave. The net change is +1 (sits) - 1 (leaves) = 0.")

    print("\nConclusion: To maximize the number of people, every move must increase the count.")
    print("This is only possible in Scenario A. Therefore, we must always place a person in a chair with no occupied neighbors.")
    print("This leads to a final state where no two occupied chairs are adjacent.")

    print("\n--- Calculation ---")
    print("The problem is now to find the maximum number of non-adjacent chairs in a row of 20.")
    print("The optimal arrangement is to fill every other chair (e.g., O E O E O...).")
    
    # For N chairs, the maximum number of non-adjacent selections is ceil(N / 2).
    max_occupied = math.ceil(num_chairs / 2)

    print(f"\nFor N = {num_chairs} chairs, the maximum number of people is calculated as:")
    # Final equation as requested
    print(f"Maximum Occupancy = ceil({num_chairs} / 2)")
    print(f"Result: {max_occupied}")

    print("\nThe numbers in the final equation are:")
    print(f"Total chairs: {num_chairs}")
    print(f"Divisor: 2")
    print(f"Maximum occupied: {max_occupied}")

    print("\nAn example of a maximum occupancy configuration is:")
    final_config = ['O' if (i % 2 == 0) else 'E' for i in range(num_chairs)]
    print("".join(final_config))

solve_chair_puzzle()
<<<10>>>