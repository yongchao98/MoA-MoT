import math

def solve_goldberg_polyhedron():
    """
    This script calculates the properties of the Goldberg polyhedron shown in the image.
    """
    
    # Step 1: Define the parameters m and n based on visual inspection.
    # The path between adjacent pentagons is m steps, a turn, and then n steps.
    # From the image, we can determine this path to be (5, 4).
    # The problem states m > n.
    m = 5
    n = 4
    
    # Step 2: The number of pentagonal faces in a Goldberg polyhedron is always 12.
    P = 12
    
    # Step 3: Calculate the number of hexagonal faces (H).
    # The formula is H = 10 * (T - 1), where T = m^2 + m*n + n^2.
    
    print(f"The identified parameters are m = {m} and n = {n}.")
    print(f"The number of pentagonal faces is always P = {P}.")
    print("\nCalculating the number of hexagonal faces (H):")
    print("The formula is H = 10 * (m^2 + m*n + n^2 - 1)")
    
    # Calculate T, the triangulation number
    t_val = m**2 + m*n + n**2
    
    # Calculate H
    h_val = 10 * (t_val - 1)
    
    # Print the equation with the numbers plugged in
    print(f"H = 10 * ({m}^2 + {m}*{n} + {n}^2 - 1)")
    print(f"H = 10 * ({m**2} + {m*n} + {n**2} - 1)")
    print(f"H = 10 * ({t_val} - 1)")
    print(f"H = 10 * ({t_val - 1})")
    print(f"H = {h_val}")
    
    # Step 4: Format and print the final answer as m,n,H,P
    print("\nThe final answer in the format m,n,H,P is:")
    print(f"{m},{n},{h_val},{P}")

solve_goldberg_polyhedron()