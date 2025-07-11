import numpy as np
from sympy.physics.wigner import wigner_3j
from sympy import N

def main():
    """
    This script solves the problem by following these steps:
    1. It recognizes that the nine images visualize Wigner 3-j symbols of the form (j j j | m1 m2 m3)
       for j = 1, 2, ..., 9.
    2. For each j, it calculates the infinity-norm of the corresponding symbol matrix. The infinity-norm
       is the maximum absolute row sum: max_{m1} sum_{m2} |(j j j | m1 m2 m3)|.
    3. It uses the `sympy` library to get the numerical values of the 3-j symbols.
    4. It iterates through j from 1 to 9, computing the norm for each case and storing it.
    5. After computing all nine norms, it finds the maximum and minimum values among them.
    6. Finally, it calculates and prints the ratio of the maximum norm to the minimum norm.
    """
    
    norms = []
    # Loop through j from 1 to 9, corresponding to the 9 images
    for j in range(1, 10):
        max_row_sum = 0.0
        # Iterate through rows, indexed by m1
        for m1 in range(-j, j + 1):
            current_row_sum = 0.0
            # Sum absolute values across the row, indexed by m2
            for m2 in range(-j, j + 1):
                m3 = -m1 - m2
                # The symbol is non-zero only if m values are valid (|m_i| <= j)
                if abs(m3) <= j:
                    # Calculate the Wigner 3-j symbol
                    symbol_val = wigner_3j(j, j, j, m1, m2, m3)
                    # Convert the sympy expression to a float and add its absolute value to the sum
                    current_row_sum += abs(float(N(symbol_val)))
            
            # Update the maximum row sum found so far for the current j
            if current_row_sum > max_row_sum:
                max_row_sum = current_row_sum
        
        norms.append(max_row_sum)

    # Find the maximum and minimum norms from the list
    max_norm = max(norms)
    min_norm = min(norms)

    # The Wigner 3-j symbol with the largest infinity-norm corresponds to j=9
    # The Wigner 3-j symbol with the smallest infinity-norm corresponds to j=2

    # Calculate the ratio
    ratio = max_norm / min_norm
    
    # Print the final result in the required format
    print("The ratio of the maximum to the minimum infinity-norm is:")
    print(f"{max_norm} / {min_norm} = {ratio}")

if __name__ == "__main__":
    # To run this code, you might need to install sympy: pip install sympy
    main()