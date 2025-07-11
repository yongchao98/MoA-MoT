import sympy
from sympy.physics.wigner import wigner_3j

def solve_wigner_ratio():
    """
    This function calculates the ratio of the maximum to the minimum infinity-norm
    among the nine Wigner 3-j symbols visualized in the image.

    The symbols are identified as (j j l; m m' mu) for j=8 and l=0, 1, ..., 8.
    The infinity-norm for a given l is max |(8 8 l; m m' -m-m')| over all m, m'.
    """

    j = 8
    norms = []

    print("Calculating the infinity-norm for each Wigner 3-j symbol (l=0 to 8)...")

    # Loop for l = 0, 1, ..., 8, which correspond to the 9 images
    for l in range(j + 1):
        max_abs_val_for_l = 0.0
        
        # Iterate over all possible magnetic quantum numbers m and m'
        for m in range(-j, j + 1):
            for mp in range(-j, j + 1):
                # The third magnetic quantum number is constrained by m1 + m2 + m3 = 0
                mu = -m - mp
                
                # Check a key selection rule: |m3| <= j3.
                # This avoids calculating symbols that are guaranteed to be zero.
                if abs(mu) <= l:
                    # Calculate the Wigner 3-j symbol using sympy
                    val = wigner_3j(j, j, l, m, mp, mu)
                    
                    # The symbol is non-zero only if it satisfies all selection rules.
                    if val != 0:
                        # Convert the symbolic result to a float and find its absolute value
                        abs_val = abs(float(val.evalf()))
                        if abs_val > max_abs_val_for_l:
                            max_abs_val_for_l = abs_val
        
        norms.append(max_abs_val_for_l)
        print(f"Infinity-norm for l={l} (Image {l+1}) is: {max_abs_val_for_l:.6f}")

    # Find the maximum and minimum norms from the list of 9 norms.
    # Filter out any potential zero norms to avoid division by zero.
    non_zero_norms = [n for n in norms if n > 1e-9] 
    
    max_total_norm = max(non_zero_norms)
    min_total_norm = min(non_zero_norms)

    # Calculate the final ratio
    ratio = max_total_norm / min_total_norm

    print("\n--- Final Calculation ---")
    print(f"Maximum infinity-norm: {max_total_norm}")
    print(f"Minimum infinity-norm: {min_total_norm}")
    print(f"Ratio = {max_total_norm} / {min_total_norm}")
    print(f"Result = {ratio}")

if __name__ == '__main__':
    solve_wigner_ratio()
