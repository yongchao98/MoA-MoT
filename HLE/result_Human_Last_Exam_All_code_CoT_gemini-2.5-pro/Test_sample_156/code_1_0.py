import sys

def solve():
    """
    Solves the alien colony puzzle by following the described plan.
    """
    # Step 1: Define initial parameters
    initial_squares = 8

    # Step 2: Explain the logic of colony growth
    print("The final size of the colony, K, is the sum of the initial 8 squares and the number of squares in any enclosed 'holes'.")
    print("A 'hole' is a vacant region completely surrounded by the initial 8 squares.")
    print(f"The {initial_squares} initial squares must form the boundary of any hole that gets filled.")
    
    # Step 3: Frame the problem as finding the maximum hole size
    print("\nThis means the number of squares on the boundary of a hole cannot be more than 8.")
    print("The problem is to find the maximum area of a hole that can be enclosed by a boundary of at most 8 squares.")

    # Step 4: Determine the maximum hole size
    # By inspection of small polyominoes (shapes on a grid):
    # - A boundary of 4 squares can enclose 1 square.
    # - A boundary of 6 squares can enclose 2 squares.
    # - A boundary of 8 squares can enclose a 2x2 area, which is 4 squares.
    # It is not possible to enclose a larger area with only 8 boundary squares.
    max_hole_area = 4
    
    print(f"\nBy inspecting possible shapes, the maximum area that can be enclosed by 8 squares is {max_hole_area}.")
    print("This can be achieved by arranging the 8 squares like this, forming a boundary around a 2x2 hole:")
    print(".XX.")
    print("X..X")
    print("X..X")
    print(".XX.")
    
    # Step 5: Calculate the final maximum size K
    max_K = initial_squares + max_hole_area
    
    print("\nStep 5: Calculate the maximal colony size, K.")
    print(f"The maximal size K is the sum of the initial squares and the maximal hole area.")
    # The final equation as requested
    print(f"K = {initial_squares} + {max_hole_area}")
    print(f"K = {max_K}")

    # Step 6: Verify that this is achievable with the given constraints
    print("\nThis configuration is achievable. For example, if d5 and e5 are the top two squares of the boundary,")
    print("the other 6 squares can be placed to form the rest of the boundary, enclosing a 2x2 hole which then gets filled.")
    print("Therefore, the maximal size of the colony is 12.")

solve()

# The final answer is wrapped in <<<>>>
sys.stdout.write("\n<<<12>>>\n")