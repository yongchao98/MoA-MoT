import math

def find_minimum_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    in a scalar field theory by checking common interaction types.
    """
    loops = 2
    min_vertices = float('inf')
    best_config = {}

    # Check for common interaction types (n-point vertices)
    # n=3 for phi^3 theory, n=4 for phi^4 theory
    for n in [3, 4]:
        # To find the minimum V, we should start with the minimum E.
        # The simplest diagrams are vacuum diagrams (E=0).
        # We check a few small even values for E.
        for E in [0, 2, 4]:
            # From L = I - V + 1 and n*V = 2*I + E, we derive:
            # I = (n*V - E)/2
            # L = (n*V - E)/2 - V + 1
            # L - 1 + E/2 = V * (n/2 - 1)
            # V = (L - 1 + E/2) / (n/2 - 1)
            
            # Avoid division by zero for n<=2
            if (n / 2 - 1) == 0:
                continue

            # Calculate the required number of vertices
            V = (loops - 1 + E / 2) / (n / 2 - 1)

            # A valid diagram must have a positive integer number of vertices.
            # We check if V is a positive integer (within a small tolerance for float issues).
            if V > 0 and abs(V - round(V)) < 1e-9:
                V = int(round(V))
                
                # Check if the corresponding number of internal lines is a non-negative integer.
                # I = (n*V - E) / 2
                I_numerator = n * V - E
                if I_numerator >= 0 and I_numerator % 2 == 0:
                    I = I_numerator // 2
                    
                    # We have found a valid configuration. Check if it's the minimum.
                    if V < min_vertices:
                        min_vertices = V
                        best_config = {'L': loops, 'V': V, 'I': I, 'E': E, 'n': n}

    # Print the results for the best configuration found
    if min_vertices != float('inf'):
        L = best_config['L']
        V = best_config['V']
        E = best_config['E']
        n = best_config['n']
        
        print("The minimum number of vertices is found using the following configuration:")
        print(f"- Interaction type: phi^{n} theory (n={n})")
        print(f"- Number of loops: L = {L}")
        print(f"- Number of external lines: E = {E} (a vacuum diagram)")
        print("\nUsing the formula: V = (L + E/2 - 1) / (n/2 - 1)")
        print("Substituting the values:")
        print(f"V = ({L} + {E}/2 - 1) / ({n}/2 - 1)")
        print(f"V = ({L - 1 + E / 2}) / ({n / 2 - 1})")
        print(f"V = {V}")
        print(f"\nTherefore, the minimum number of vertices is {V}.")

    else:
        print("No valid configuration found.")

find_minimum_vertices()
<<<1>>>