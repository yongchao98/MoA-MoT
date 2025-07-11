import numpy as np
from collections import Counter

def find_coordination_polyhedra():
    """
    This script parses CIF-like crystal data, finds the nearest neighbor
    coordination for each unique atom, determines the polyhedron formula,
    and prints the result in the specified format.
    """
    # 1. Parse CIF data for ReAl12
    cell_lengths = {'a': 7.609100, 'b': 6.611700, 'c': 9.023000}
    
    unique_sites = [
        {'label': 'Al_A', 'symbol': 'Al', 'fract': np.array([0.3182, 0.2158, 0.2500])},
        {'label': 'Al_B', 'symbol': 'Al', 'fract': np.array([0.0000, 0.3662, 0.1030])},
        {'label': 'Al_C', 'symbol': 'Al', 'fract': np.array([0.1743, 0.0000, 0.0000])},
        {'label': 'Re_A', 'symbol': 'Re', 'fract': np.array([0.0000, 0.0445, 0.2500])}
    ]

    sym_ops_str = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2',
        '-x+1/2,y+1/2,z', '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2',
        'x+1/2,y+1/2,z', 'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
    ]

    # 2. Utility functions
    def parse_symop(op_str):
        """Parses a symmetry operation string into a callable function."""
        parts = op_str.split(',')
        # Using eval is concise for this controlled input.
        return lambda x, y, z: (eval(parts[0]), eval(parts[1]), eval(parts[2]))

    sym_ops = [parse_symop(s) for s in sym_ops_str]

    a, b, c = cell_lengths['a'], cell_lengths['b'], cell_lengths['c']
    def frac_to_cart(fract_coords):
        """Converts fractional coordinates to Cartesian for an orthorhombic cell."""
        return np.array([fract_coords[0] * a, fract_coords[1] * b, fract_coords[2] * c])

    # 3. & 4. Find coordination polyhedra for each unique site
    polyhedra_results = []
    for center_site in unique_sites:
        center_label_base = center_site['symbol']
        center_fract = center_site['fract']
        
        distances = []
        
        # Find all neighbors using the Minimum Image Convention
        for neighbor_site in unique_sites:
            for op in sym_ops:
                sym_neighbor_fract = np.array(op(neighbor_site['fract'][0], neighbor_site['fract'][1], neighbor_site['fract'][2]))
                
                # Vector from center to neighbor in fractional coordinates
                diff = sym_neighbor_fract - center_fract
                # Apply periodic boundary conditions (Minimum Image Convention)
                diff -= np.round(diff)
                
                dist = np.linalg.norm(frac_to_cart(diff))
                
                # Exclude the central atom itself (distance ~ 0)
                if dist > 1e-4:
                    distances.append({
                        'dist': dist,
                        'symbol': neighbor_site['symbol']
                    })

        # Sort neighbors by distance to find shells
        distances.sort(key=lambda x: x['dist'])

        # Identify the first coordination shell by finding a significant distance gap
        if not distances:
            continue
            
        shell_cutoff_index = len(distances) # Default to all if no gap is found
        for i in range(len(distances) - 1):
            # A gap is defined as a distance increase of more than 20%
            if distances[i+1]['dist'] / distances[i]['dist'] > 1.20:
                shell_cutoff_index = i + 1
                break
        
        first_shell = distances[:shell_cutoff_index]

        # 5. Determine the formula of the polyhedron
        shell_composition = Counter(d['symbol'] for d in first_shell)
        
        formula_parts = []
        for symbol in sorted(shell_composition.keys()):
            count = shell_composition[symbol]
            formula_parts.append(f"{symbol}{count if count > 1 else ''}")
        
        formula = "".join(formula_parts)
        polyhedra_results.append(f"{center_label_base}, {formula}")

    # 6. Format the final output, removing duplicate environments
    unique_results = sorted(list(set(polyhedra_results)))
    print("; ".join(unique_results))

find_coordination_polyhedra()