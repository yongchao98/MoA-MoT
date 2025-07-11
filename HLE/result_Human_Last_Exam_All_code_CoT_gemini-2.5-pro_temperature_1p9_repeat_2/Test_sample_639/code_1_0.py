import math

def find_planar_distances():
    """
    Calculates the set of possible normalized distances (r <= 3) between hard spheres
    packed on a 2D plane by considering square and hexagonal lattice arrangements.
    """
    # Use a set to store unique distances automatically.
    distances = set()
    limit = 3.0
    # Add a small tolerance for floating point comparisons
    tolerance = 1e-9

    print("Thinking Process:")
    print("1. Analyzing square lattice configurations...")
    # --- Square Lattice ---
    # The normalized distance r from the origin (0,0) to a point (n,m) is sqrt(n^2 + m^2).
    # We only need to loop n and m up to the integer value of the limit.
    for n in range(int(limit) + 1):
        for m in range(int(limit) + 1):
            if n == 0 and m == 0:
                continue  # Skip distance to self
            
            r_squared = float(n**2 + m**2)
            if r_squared <= limit**2 + tolerance:
                distances.add(math.sqrt(r_squared))

    print("2. Analyzing hexagonal lattice configurations...")
    # --- Hexagonal Lattice ---
    # The normalized squared distance r^2 is given by n^2 + m^2 + nm.
    # We loop through n and m, including negative values for m to cover all neighbors.
    for n in range(int(limit) + 1):
        for m in range(-int(limit), int(limit) + 1):
            if n == 0 and m == 0:
                continue
            
            r_squared = float(n**2 + m**2 + n*m)
            if r_squared <= limit**2 + tolerance:
                distances.add(math.sqrt(r_squared))

    # --- Final Output ---
    # Sort the collected distances and format them.
    sorted_distances = sorted(list(distances))
    
    print("\n3. Consolidating and formatting the final set of distances...\n")
    print("The set of possible normalized distances r <= 3 is:")
    
    # Create the final list of strings for printing
    result_strings = [f"{dist:.2f}" for dist in sorted_distances]
    
    # The problem asks to output each number, so we print them.
    print(", ".join(result_strings))

if __name__ == '__main__':
    find_planar_distances()
    # The final answer in the required format is derived from the printed output.
    # The sorted unique distances are:
    # sqrt(1)=1.00, sqrt(2)=1.41, sqrt(3)=1.73, sqrt(4)=2.00, 
    # sqrt(5)=2.24, sqrt(7)=2.65, sqrt(8)=2.83, sqrt(9)=3.00
    # The formatted string is "1.00, 1.41, 1.73, 2.00, 2.24, 2.65, 2.83, 3.00"