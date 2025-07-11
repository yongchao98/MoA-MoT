import sys

def find_minimum_vertices():
    """
    This function searches for the minimum number of vertices in a two-loop
    Feynman diagram for a generic interacting scalar field theory.
    
    The number of vertices (V) is related to the interaction order (n) and
    the number of external lines (E) by the formula:
    V = (2 + E) / (n - 2)
    
    We search for the minimum integer V by checking different values for n and E.
    Constraints:
    - n >= 3 (for an interaction)
    - E >= 0
    - V must be a positive integer
    """
    
    min_V = float('inf')
    best_n = -1
    best_E = -1

    # Search over a plausible range of n and E
    # We prioritize lower E and check increasing n
    for E_val in range(10): # Number of external lines
        for n_val in range(3, 10): # Order of the interaction vertex
            
            numerator = 2 + E_val
            denominator = n_val - 2
            
            # V must be an integer, so the numerator must be divisible by the denominator
            if numerator % denominator == 0:
                V = numerator // denominator
                
                # We are looking for the minimum positive V
                if 0 < V < min_V:
                    min_V = V
                    best_n = n_val
                    best_E = E_val

    # Print the result and the calculation for the optimal case
    print("The minimum number of vertices is found using the general formula: V = (2 + E) / (n - 2)")
    print(f"The search has found the minimum for a phi^{best_n} theory (n={best_n}) and a diagram with {best_E} external lines (E={best_E}).")
    print("\nFinal Calculation:")
    print(f"V = (2 + {best_E}) / ({best_n} - 2)")
    print(f"V = ({2 + best_E}) / ({best_n - 2})")
    print(f"V = {min_V}")
    print("\nSo, the minimum number of vertices in a two-loop Feynman diagram for an interacting scalar field theory is {min_V}.")
    

# Execute the function
find_minimum_vertices()