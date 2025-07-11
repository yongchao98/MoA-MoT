import sympy.physics.wigner
import sys

def calculate_infinity_norm(j1, j2, j3):
    """
    Calculates the infinity-norm (maximum absolute value) of a Wigner 3-j symbol.
    The symbol is calculated for all valid m1 and m2.
    """
    max_abs_val = 0.0
    # m-values range from -j to +j
    m1_range = range(-j1, j1 + 1)
    m2_range = range(-j2, j2 + 1)

    for m1 in m1_range:
        for m2 in m2_range:
            # The third m-quantum number is fixed by the other two
            m3 = -m1 - m2
            
            # The symbol is non-zero only if |m3| <= j3
            if abs(m3) <= j3:
                # Calculate the 3-j symbol using sympy
                val = sympy.physics.wigner.wigner_3j(j1, j2, j3, m1, m2, m3)
                
                # Update the maximum absolute value found so far
                if val != 0:
                    # .evalf() evaluates the symbolic expression to a float
                    current_abs_val = abs(float(val.evalf()))
                    if current_abs_val > max_abs_val:
                        max_abs_val = current_abs_val
    return max_abs_val

def solve_wigner_ratio():
    """
    Main function to solve the problem.
    """
    print("Identifying parameters for the 9 Wigner 3-j symbols...")
    # These parameters correspond to the 9 plots in the image.
    j1_values = range(8, 17)  # j1 = 8, 9, 10, ..., 16
    j2 = 8
    j3 = 8
    print(f"Parameters identified: j2={j2}, j3={j3}, and j1 varies from 8 to 16.\n")
    
    norms = []
    print("Calculating the infinity-norm for each symbol... (This may take a moment)")
    
    # Loop through each of the 9 symbols
    for j1 in j1_values:
        # The triangle inequality for j-values (j1, j2, j3) must be satisfied.
        # |j2 - j3| <= j1 <= j2 + j3  =>  |8-8| <= j1 <= 8+8  =>  0 <= j1 <= 16
        # All our j1 values (8-16) satisfy this.
        
        sys.stdout.write(f"Calculating for (j1={j1}, j2={j2}, j3={j3})... ")
        sys.stdout.flush()
        
        norm = calculate_infinity_norm(j1, j2, j3)
        norms.append(norm)
        
        print(f"Norm = {norm:.6f}")

    # Find the maximum and minimum norms from the list
    max_norm = max(norms)
    min_norm = min(norms)
    
    # Calculate the final ratio
    if min_norm == 0:
        ratio = float('inf')
    else:
        ratio = max_norm / min_norm
        
    print("\n--- Calculation Complete ---")
    print(f"The nine calculated infinity-norms are:")
    for i, norm in enumerate(norms):
        print(f"  Symbol with j1={j1_values[i]}: {norm:.6f}")
        
    print("\n--- Final Answer ---")
    print(f"Maximum norm found: {max_norm}")
    print(f"Minimum norm found: {min_norm}")
    print(f"The ratio of maximum to minimum norm is:")
    print(f"{max_norm} / {min_norm} = {ratio}")
    
    return ratio

if __name__ == '__main__':
    final_ratio = solve_wigner_ratio()
    # The final answer is wrapped in <<<>>>
    print(f"\n<<<{final_ratio}>>>")