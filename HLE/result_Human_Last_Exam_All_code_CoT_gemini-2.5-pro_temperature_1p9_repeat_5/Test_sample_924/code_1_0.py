def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes
    in a 3x3x3 cube with the given face constraints.
    """
    print("Finding the smallest and largest number of green cubes.")
    print("-" * 50)
    
    min_total_green = float('inf')
    max_total_green = float('-inf')

    # Gf is the number of green face-center cubes. It must be a multiple of 3.
    possible_Gf = [0, 3, 6]
    
    # G_core can be 0 (red) or 1 (green).
    possible_G_core = [0, 1]

    print("Step 1: Derived relations between numbers of green cubes by type:")
    print("3*Gc + 2*Ge + Gf = 36")
    print("Ge + Gf = 12")
    print("From these, we get: 3*Gc - Gf = 12\n")

    for Gf in possible_Gf:
        # From 3*Gc - Gf = 12  =>  Gc = (12 + Gf) / 3
        Gc = (12 + Gf) / 3
        if Gc.is_integer():
            Gc = int(Gc)
            # From Ge + Gf = 12 => Ge = 12 - Gf
            Ge = 12 - Gf
            
            print(f"Case: Number of green face-centers (Gf) = {Gf}")
            print(f"  - Number of green corners (Gc) = (12 + {Gf}) / 3 = {Gc}")
            print(f"  - Number of green edges (Ge) = 12 - {Gf} = {Ge}")

            for G_core in possible_G_core:
                G_total = Gc + Ge + Gf + G_core
                
                print(f"    - Core cube is {'green' if G_core == 1 else 'red'} (G_core = {G_core})")
                print(f"    - Total Green = Gc + Ge + Gf + G_core = {Gc} + {Ge} + {Gf} + {G_core} = {G_total}")
                
                if G_total < min_total_green:
                    min_total_green = G_total
                if G_total > max_total_green:
                    max_total_green = G_total
            print("-" * 20)
            
    print("\n" + "=" * 50)
    print("Summary of Results:")
    print(f"The smallest possible number of green cubes is: {min_total_green}")
    print(f"The largest possible number of green cubes is: {max_total_green}")
    print("=" * 50)

solve_cube_problem()
<<<16, 19>>>