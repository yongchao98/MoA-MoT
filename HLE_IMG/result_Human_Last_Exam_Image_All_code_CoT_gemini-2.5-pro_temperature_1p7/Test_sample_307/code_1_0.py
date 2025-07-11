import sympy
from sympy.physics.wigner import wigner_3j

def solve_wigner_ratio():
    """
    Calculates the ratio of the maximum to the minimum infinity-norm for the
    nine Wigner 3-j symbols visualized in the image.
    """
    # The parameters for the 3-j symbols are (j1=8, j2=8, j3)
    j1 = 8
    j2 = 8
    j3_values = range(8, 17)  # Corresponds to plots 1 through 9 for j3=8..16

    infinity_norms = []
    print("Calculating the infinity-norm for each of the 9 plots...")

    # Calculate the infinity norm for each j3
    for j3 in j3_values:
        max_abs_val = 0.0
        # Iterate over all allowed magnetic quantum numbers
        for m1 in range(-j1, j1 + 1):
            for m3 in range(-j3, j3 + 1):
                # The third magnetic quantum number is fixed by m1+m2+m3=0
                m2 = -m1 - m3
                
                # The symbol is zero unless |m2| <= j2
                if abs(m2) <= j2:
                    # Use sympy to calculate the Wigner 3-j symbol
                    val = wigner_3j(j1, j2, j3, m1, m2, m3)
                    # The result is symbolic, evaluate it to a float
                    val_float = float(val.evalf())
                    
                    # Update the maximum absolute value for the current j3
                    if abs(val_float) > max_abs_val:
                        max_abs_val = abs(val_float)
        
        infinity_norms.append(max_abs_val)
        plot_number = j3 - 7
        print(f"Plot {plot_number} (j3={j3}): Infinity-norm = {max_abs_val:.8f}")

    # Find the maximum and minimum norms from the list
    max_norm = max(infinity_norms)
    min_norm = min(infinity_norms)

    # Calculate the final ratio
    ratio = max_norm / min_norm
    
    print("\n--- Final Result ---")
    print(f"The maximum infinity-norm found is: {max_norm:.8f}")
    print(f"The minimum infinity-norm found is: {min_norm:.8f}")
    print(f"\nThe ratio of the maximum to the minimum infinity-norm is:")
    print(f"{max_norm:.8f} / {min_norm:.8f} = {ratio:.8f}")

if __name__ == '__main__':
    solve_wigner_ratio()