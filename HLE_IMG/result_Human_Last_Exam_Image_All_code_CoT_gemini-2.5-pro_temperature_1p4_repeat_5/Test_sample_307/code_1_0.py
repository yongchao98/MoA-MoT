import numpy as np
from sympy.physics.wigner import wigner_3j
from sympy import N

def calculate_wigner_infinity_norm(j1, j2, j3):
    """
    Calculates the infinity-norm for a matrix of Wigner 3-j symbols
    for a given (j1, j2, j3). The matrix is indexed by (m1, m2).
    A 'row' corresponds to a fixed m2.
    """
    max_abs_row_sum = 0
    
    # Iterate over each row (each value of m2)
    for m2 in range(-j2, j2 + 1):
        current_abs_row_sum = 0
        # Sum the absolute values of the symbols across the row (all m1 values)
        for m1 in range(-j1, j1 + 1):
            m3 = -m1 - m2
            
            # The symbol is non-zero only if |m3| <= j3.
            # This is a fundamental selection rule.
            if abs(m3) <= j3:
                # Use SymPy's wigner_3j function and convert the result to a float.
                symbol_value = N(wigner_3j(j1, j2, j3, m1, m2, m3))
                current_abs_row_sum += abs(symbol_value)
                
        # Keep track of the maximum sum found so far.
        if current_abs_row_sum > max_abs_row_sum:
            max_abs_row_sum = current_abs_row_sum
            
    return float(max_abs_row_sum)

def main():
    """
    Main function to calculate and print the required ratio.
    """
    # Parameters identified from the problem description and image source.
    j3 = 8
    j_pairs = [
        (7, 7), (7, 8), (7, 9),
        (8, 7), (8, 8), (8, 9),
        (9, 7), (9, 8), (9, 9)
    ]
    
    print("Calculating infinity-norm for each of the 9 matrices...\n")
    
    infinity_norms = []
    for i, (j1, j2) in enumerate(j_pairs):
        norm = calculate_wigner_infinity_norm(j1, j2, j3)
        infinity_norms.append(norm)
        print(f"Norm for matrix {i+1} (j1={j1}, j2={j2}, j3={j3}): {norm}")
        
    # Find the maximum and minimum of the calculated norms.
    max_norm = max(infinity_norms)
    min_norm = min(infinity_norms)
    
    # Calculate the final ratio.
    ratio = max_norm / min_norm
    
    print("\n--- Summary ---")
    print(f"Maximum infinity-norm: {max_norm}")
    print(f"Minimum infinity-norm: {min_norm}")
    
    print("\nFinal Equation:")
    print(f"{max_norm} / {min_norm} = {ratio}")

if __name__ == "__main__":
    main()