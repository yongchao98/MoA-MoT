import math

def calculate_and_print_distances():
    """
    Calculates and prints the possible normalized distances (r) between hard spheres
    in a 2D jammed configuration for r <= 3.

    The calculation is based on a triangular lattice model, which represents the
    densest and most stable packing of disks on a plane. The position of any
    sphere relative to a central sphere at the origin can be described by
    integer linear combinations (m, n) of two basis vectors. The squared distance
    is given by the formula: r^2 = m^2 + m*n + n^2.
    """
    print("Calculating possible distances r <= 3 for 2D packed hard spheres.")
    print("The formula for squared distance is r^2 = m^2 + m*n + n^2")
    print("-" * 60)

    # We search for integer pairs (m, n) where r^2 <= 9 and find the unique distances.
    # The resulting r^2 values are 1, 3, 4, 7, 9.
    
    # We choose one representative (m, n) pair for each unique distance.
    configurations = [
        {'m': 1, 'n': 0, 'desc': 'Two spheres in direct contact.'},
        {'m': 1, 'n': 1, 'desc': 'Two spheres separated by one sphere in a triangular well.'},
        {'m': 2, 'n': 0, 'desc': 'Two spheres separated by one sphere in a straight line.'},
        {'m': 2, 'n': 1, 'desc': 'A next-nearest neighbor configuration in the lattice.'},
        {'m': 3, 'n': 0, 'desc': 'Two spheres separated by two spheres in a straight line.'}
    ]

    final_distances = []

    for config in configurations:
        m, n = config['m'], config['n']
        r_squared = m**2 + m*n + n**2
        r = math.sqrt(r_squared)
        final_distances.append(r)
        
        print(f"Distance from configuration (m={m}, n={n}):")
        print(f"   Description: {config['desc']}")
        print(f"   Equation: r = sqrt({m}^2 + {m}*{n} + {n}^2) = sqrt({r_squared})")
        print(f"   Result: r = {r:.2f}\n")

    print("-" * 60)
    print("The set of unique distances r <= 3 is:")
    
    # Sort and format the final list for clear output
    sorted_distances = sorted(list(set(final_distances)))
    formatted_distances = [f"{dist:.2f}" for dist in sorted_distances]
    
    print(formatted_distances)

calculate_and_print_distances()

# The final answer is the set of calculated distances
answer = sorted(list(set(d for d in [math.sqrt(cfg['m']**2 + cfg['m']*cfg['n'] + cfg['n']**2) for cfg in [
        {'m': 1, 'n': 0}, {'m': 1, 'n': 1}, {'m': 2, 'n': 0}, {'m': 2, 'n': 1}, {'m': 3, 'n': 0}
    ]])))
formatted_answer = ', '.join([f'{x:.2f}' for x in answer])
print(f'<<<{formatted_answer}>>>')