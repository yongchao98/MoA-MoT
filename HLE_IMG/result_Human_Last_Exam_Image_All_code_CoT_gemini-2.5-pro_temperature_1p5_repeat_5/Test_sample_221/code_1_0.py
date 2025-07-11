import math

def solve_goldberg():
    """
    This function determines the parameters of the Goldberg polyhedron from the image and calculates the number of faces.
    
    The steps are:
    1. Identify the number of pentagonal faces, P, which is always 12 for a Goldberg polyhedron.
    2. Visually determine the parameters m and n by tracing the path between adjacent pentagons.
       The path involves m steps, a 60-degree turn, and n steps. From the image, m=5 and n=3.
    3. Calculate the number of hexagonal faces, H, using the formula H = 10 * (m^2 + m*n + n^2 - 1).
    4. Print the final answer in the format m,n,H,P.
    """
    
    # Parameters m and n are determined by visual inspection of the path between pentagons.
    # The path consists of 5 steps, a turn, and 3 steps. As per the hint m > n.
    m = 5
    n = 3
    
    # The number of pentagonal faces (P) in any Goldberg polyhedron is constant.
    P = 12
    
    # Calculate the number of hexagonal faces (H) using the standard formula.
    H = 10 * (m**2 + m * n + n**2 - 1)
    
    # Output the result in the specified format: m,n,H,P
    print(f"{m},{n},{int(H)},{P}")

solve_goldberg()