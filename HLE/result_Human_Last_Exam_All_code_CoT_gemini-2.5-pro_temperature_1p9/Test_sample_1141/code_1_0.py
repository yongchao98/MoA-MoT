import math

def solve_cutting_problem():
    """
    Analyzes the cutting stock problem to find the optimal solution.
    """
    # --- Problem Definition ---
    billet_cm = (16, 11, 4)
    prices = {'B2': 150, 'T1': 5, 'B1': 1}
    dims = {
        'B2': {'radius_cm': 2},
        'T1': {'side_cm': 1},
        'B1': {'diameter_cm': 1, 'radius_cm': 0.5}
    }

    # --- Analysis ---
    print("Step 1: Analyzing Value per Unit Volume")
    
    # Calculate volumes in cm^3
    vol_b2 = (4/3) * math.pi * (dims['B2']['radius_cm'] ** 3)
    vol_t1 = dims['T1']['side_cm'] ** 3
    vol_b1 = (4/3) * math.pi * (dims['B1']['radius_cm'] ** 3)
    
    # Calculate value density (value per cm^3)
    val_density_b2 = prices['B2'] / vol_b2
    val_density_t1 = prices['T1'] / vol_t1
    val_density_b1 = prices['B1'] / vol_b1
    
    print(f"Value density of B2 (2cm radius ball): {val_density_b2:.2f} per cm^3")
    print(f"Value density of T1 (1cm cube): {val_density_t1:.2f} per cm^3")
    print(f"Value density of B1 (1cm diameter ball): {val_density_b1:.2f} per cm^3")
    print("\nObservation: T1 cubes offer the highest value per unit of volume.")
    print("Because cubes can be packed perfectly (tile space), this suggests a T1-only strategy might be optimal.\n")

    print("Step 2: Evaluating a T1-only Strategy")
    
    num_t1_x = int(billet_cm[0] / dims['T1']['side_cm'])
    num_t1_y = int(billet_cm[1] / dims['T1']['side_cm'])
    num_t1_z = int(billet_cm[2] / dims['T1']['side_cm'])
    total_t1_possible = num_t1_x * num_t1_y * num_t1_z
    
    t1_only_value = total_t1_possible * prices['T1']
    
    print(f"The 16x11x4 cm billet can be perfectly packed with {num_t1_x}x{num_t1_y}x{num_t1_z} = {total_t1_possible} T1 cubes.")
    print(f"Total value from this strategy: {total_t1_possible} * {prices['T1']} = {t1_only_value}\n")

    print("Step 3: Checking if Replacing T1s with B2s is Profitable")
    # A B2 ball has a 4cm diameter. Its bounding box is 4x4x4 cm.
    # To place one B2, we must remove the T1s occupying at least this volume.
    t1s_removed_for_b2 = (2 * dims['B2']['radius_cm'])**3 / vol_t1
    value_lost_from_t1 = t1s_removed_for_b2 * prices['T1']
    value_gained_from_b2 = prices['B2']
    net_change = value_gained_from_b2 - value_lost_from_t1
    
    print("To place one B2 ball (diameter 4cm), we must clear a space of at least 4x4x4 cm.")
    print(f"This removes {int(t1s_removed_for_b2)} T1 cubes, a value loss of {int(value_lost_from_t1)}.")
    print(f"The B2 ball only adds a value of {value_gained_from_b2}.")
    print(f"The net change in value is {int(net_change)}, which is a loss. Therefore, it is not profitable to replace T1 cubes with B2 balls.")
    print("A similar analysis shows it is also not profitable to replace T1 cubes with B1 balls.\n")
    
    print("--- Conclusion ---")
    print("The optimal strategy is to cut the billet entirely into T1 cubes, maximizing the total value.")
    
    final_counts = {'B2': 0, 'T1': total_t1_possible, 'B1': 0}
    final_total_value = t1_only_value
    
    print("\nFinal Optimal Cutting Plan:")
    print(f"Pieces of type B2: {final_counts['B2']}")
    print(f"Pieces of type T1: {final_counts['T1']}")
    print(f"Pieces of type B1: {final_counts['B1']}")
    
    print("\nFinal Equation (Count * Price):")
    print(f"({final_counts['B2']} * {prices['B2']}) + ({final_counts['T1']} * {prices['T1']}) + ({final_counts['B1']} * {prices['B1']}) = {final_total_value}")

if __name__ == '__main__':
    solve_cutting_problem()