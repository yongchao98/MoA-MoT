import math

def solve_for_min_vertices():
    """
    This script calculates the minimum number of vertices for a two-loop
    Feynman diagram in an interacting scalar field theory.
    """
    
    print("To find the minimum number of vertices in a two-loop diagram, we use two fundamental relations for Feynman diagrams.")
    print("Let V = number of vertices, I = internal lines, E = external lines, L = loops, and n = interaction order (e.g., n=4 for φ^4 theory).")
    print("\nThe two key formulas are:")
    print("1. Loop Formula:      L = I - V + 1")
    print("2. Vertex Formula:    n * V = 2 * I + E")
    print("-" * 50)
    
    # We are given the number of loops, L = 2.
    L = 2
    
    print(f"Given L = {L}, we can rearrange the Loop Formula to solve for I:")
    print(f"{L} = I - V + 1   =>   I = V + {L} - 1   =>   I = V + 1")
    print("\nNow, substitute this expression for I into the Vertex Formula:")
    print("n * V = 2 * (V + 1) + E")
    print("n * V = 2*V + 2 + E")
    print("(n - 2) * V = 2 + E")
    print("\nThis gives our final formula for the number of vertices:")
    print("V = (2 + E) / (n - 2)")
    print("-" * 50)
    
    print("To find the minimum positive integer V, we must test combinations of integer n (≥ 3) and E (≥ 0).")

    min_V = float('inf')
    min_params = None

    # We iterate through possible values of n and E to find the minimum valid V.
    for n in range(3, 10):  # Test for theories like φ^3, φ^4, ...
        for E in range(0, 10): # Test for different numbers of external particles
            
            # The denominator (n-2) must be positive, which is true for n>=3.
            numerator = L - 1 + L - 1 + E  # Generalized formula: V = (2L - 2 + E) / (n - 2)
            numerator = 2 + E              # Simplified for L=2
            denominator = n - 2
            
            # V must be a positive integer.
            if numerator % denominator == 0:
                V = numerator // denominator
                if V > 0 and V < min_V:
                    # nV and E must have the same parity since nV - E = 2I (which is even).
                    if (n * V) % 2 == E % 2:
                        min_V = V
                        min_params = {'n': n, 'E': E, 'V': V}

    if min_params:
        n_min = min_params['n']
        E_min = min_params['E']
        V_min = min_params['V']
        
        print(f"\nThe search finds a minimum when V = {V_min}.")
        print(f"This occurs for several combinations, the simplest being for a φ^{n_min} theory with E = {E_min} external lines (a vacuum diagram).")
        print("\nPlugging these values into the formula gives the final calculation:")
        
        # The prompt requires outputting each number in the final equation.
        print(f"\n{V_min} = ({L + L - 2} + {E_min}) / ({n_min} - 2)")

    else:
        print("Could not find a valid combination in the tested range.")

if __name__ == '__main__':
    solve_for_min_vertices()
    # The final answer as deduced from the logic above
    print("\n<<<1>>>")