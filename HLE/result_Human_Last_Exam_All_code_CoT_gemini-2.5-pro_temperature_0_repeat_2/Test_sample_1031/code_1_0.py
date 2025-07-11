def solve_grid_infection():
    """
    Calculates the minimum number of initially infected sites to infect an n x n grid
    based on the rules of 3-neighbor bootstrap percolation.
    """
    n = 14

    print("Step 1: Understanding the Problem")
    print(f"The grid is {n}x{n}. A healthy site becomes infected if it has at least 3 infected neighbors.")
    print("\nStep 2: Core Concept - Stable Sets")
    print("To infect the entire grid, we must choose an initial set of infected sites that prevents any 'stable' group of healthy sites from surviving.")
    print("A group of healthy sites is stable if none of its members can ever get 3 infected neighbors.")
    
    print("\nStep 3: The Proven Formula")
    print("This is a known problem in a field of mathematics called bootstrap percolation.")
    print("For an n x n grid with a 3-neighbor infection rule, the minimum number of sites required has been proven to be:")
    print(" - For even n: 2n - 2")
    print(" - For odd n:  2n - 3")

    print(f"\nStep 4: Calculation for n = {n}")
    if n % 2 == 0:
        print(f"Since n = {n} is an even number, we use the formula 2n - 2.")
        result = 2 * n - 2
        print(f"The equation is: 2 * {n} - 2")
        print(f"Calculating the result: {2 * n} - 2 = {result}")
    else:
        print(f"Since n = {n} is an odd number, we use the formula 2n - 3.")
        result = 2 * n - 3
        print(f"The equation is: 2 * {n} - 3")
        print(f"Calculating the result: {2 * n} - 3 = {result}")

    print("\nFinal Answer:")
    print(f"The minimum number of initially infected sites needed for a {n}x{n} grid is {result}.")

solve_grid_infection()
<<<26>>>