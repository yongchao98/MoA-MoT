import math

def find_2d_packing_distances():
    """
    Calculates the possible normalized distances (r <= 3) between hard spheres
    on a 2D hexagonal lattice.
    
    The squared distance r^2 is given by r^2 = n^2 + n*m + m^2, where n and
    m are integers defining the lattice position.
    """
    
    print("The possible normalized distances r between sphere centers in a 2D dense packing for r <= 3 are:")
    
    # Use a dictionary to store the result and the (n, m) pair that generated it
    # to avoid duplicates and show a representative equation.
    # Key: squared distance, Value: (n, m, formatted_result_string)
    distances_map = {}
    
    # We are looking for r <= 3, which means r^2 <= 9.
    # We can limit the search space for n and m.
    # If n=0, m^2 <= 9 -> |m|<=3
    # If n=1, 1+m+m^2 <= 9 -> m^2+m-8<=0 -> m is approx between -3.3 and 2.3
    # If n=2, 4+2m+m^2 <= 9 -> m^2+2m-5<=0 -> m is approx between -3.4 and 1.4
    # If n=3, 9+3m+m^2 <= 9 -> m^2+3m<=0 -> m is between -3 and 0
    # So, n in [0, 3] and m in [-3, 3] is a safe range to check.
    
    limit = 3
    for n in range(limit + 1):
        for m in range(-limit, limit + 1):
            
            # This combination represents the central sphere itself
            if n == 0 and m == 0:
                continue

            r_squared = n**2 + n * m + m**2
            
            if r_squared <= limit**2:
                # If we haven't seen this distance yet, store it
                if r_squared not in distances_map:
                    r = math.sqrt(r_squared)
                    # We store the result and a sample calculation
                    # For negative m, show the calculation correctly
                    m_sign = "-" if m < 0 else "+"
                    m_abs = abs(m)
                    equation = f"r = sqrt({n}^2 + {n}*({m}) + ({m})^2) = sqrt({r_squared})"
                    # Simplified for printing
                    simple_r = "sqrt({})".format(r_squared) if not r.is_integer() else int(r)
                    result_string = f"{equation} = {r:.2f}"
                    distances_map[r_squared] = (n, m, result_string)

    # Sort the distances by their value
    sorted_distances = sorted(distances_map.keys())

    # Print the final results
    for r_sq in sorted_distances:
        print(distances_map[r_sq][2])

find_2d_packing_distances()

# The final result is the set of calculated distances
final_distances = [1.00, 1.73, 2.00, 2.65, 3.00]
# The user wants the answer in the format <<<answer content>>>.
# Based on the problem, the set of real numbers is requested.
print(f"<<<{final_distances}>>>")