import math
from collections import Counter

def solve_polyhedra():
    """
    Analyzes a crystal structure from CIF data to find coordination polyhedra.
    """
    # 1. Parse Crystal Data (hardcoded from the problem description)
    cell_lengths = {'a': 7.609100, 'b': 6.611700, 'c': 9.023000}
    
    # Asymmetric unit atoms
    asym_atoms = [
        {'label': 'Al_A', 'symbol': 'Al', 'x': 0.3182, 'y': 0.2158, 'z': 0.2500},
        {'label': 'Al_B', 'symbol': 'Al', 'x': 0.0000, 'y': 0.3662, 'z': 0.1030},
        {'label': 'Al_C', 'symbol': 'Al', 'x': 0.1743, 'y': 0.0000, 'z': 0.0000},
        {'label': 'Re_A', 'symbol': 'Re', 'x': 0.0000, 'y': 0.0445, 'z': 0.2500},
    ]

    sym_ops_str = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2',
        '-x+1/2,y+1/2,z', '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2',
        'x+1/2,y+1/2,z', 'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
    ]

    def apply_sym_op(op_str, x, y, z):
        """Applies a symmetry operation to a coordinate."""
        new_x, new_y, new_z = op_str.split(',')
        # Safely evaluate the expression
        local_vars = {'x': x, 'y': y, 'z': z}
        return (eval(new_x, {}, local_vars),
                eval(new_y, {}, local_vars),
                eval(new_z, {}, local_vars))

    # 2. Generate a Full Unit Cell
    unit_cell_atoms = []
    seen_coords = set()
    for atom in asym_atoms:
        for op in sym_ops_str:
            nx, ny, nz = apply_sym_op(op, atom['x'], atom['y'], atom['z'])
            nx, ny, nz = nx % 1, ny % 1, nz % 1  # Keep coordinates in [0, 1)
            coord_tuple = (round(nx, 4), round(ny, 4), round(nz, 4))
            if coord_tuple not in seen_coords:
                unit_cell_atoms.append({'symbol': atom['symbol'], 'x': nx, 'y': ny, 'z': nz})
                seen_coords.add(coord_tuple)
    
    # 3. Construct a Supercell
    supercell_atoms = []
    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                for atom in unit_cell_atoms:
                    supercell_atoms.append({
                        'symbol': atom['symbol'],
                        'x': atom['x'] + i, 'y': atom['y'] + j, 'z': atom['z'] + k
                    })

    def cartesian_distance(p1, p2, dims):
        """Calculates distance between two points in orthorhombic cell."""
        dx = (p1['x'] - p2['x']) * dims['a']
        dy = (p1['y'] - p2['y']) * dims['b']
        dz = (p1['z'] - p2['z']) * dims['c']
        return math.sqrt(dx**2 + dy**2 + dz**2)

    unique_polyhedra = {}
    GAP_TOLERANCE = 1.25  # A >25% jump in distance indicates a new shell

    # 4 & 5. Find Nearest Neighbors and Identify Coordination Shell
    for central_atom in asym_atoms:
        distances = []
        for neighbor_atom in supercell_atoms:
            dist = cartesian_distance(central_atom, neighbor_atom, cell_lengths)
            if dist > 0.001:  # Exclude the central atom itself
                distances.append({'dist': dist, 'symbol': neighbor_atom['symbol']})
        
        distances.sort(key=lambda item: item['dist'])

        first_shell_neighbors = []
        if not distances:
            continue
        
        first_dist = distances[0]['dist']
        first_shell_neighbors.append(distances[0]['symbol'])

        for i in range(1, len(distances)):
            if distances[i]['dist'] > first_dist * GAP_TOLERANCE:
                break
            first_shell_neighbors.append(distances[i]['symbol'])

        # 6. Determine Polyhedron Formula
        polyhedron_counts = Counter(first_shell_neighbors)
        # Sort elements alphabetically for a consistent formula
        formula = "".join(
            f"{elem}{polyhedron_counts[elem]}" if polyhedron_counts[elem] > 1 else elem
            for elem in sorted(polyhedron_counts.keys())
        )
        
        unique_polyhedra[central_atom['symbol']] = formula

    # 7. Format the Output
    output_parts = []
    for symbol, formula in sorted(unique_polyhedra.items()):
        output_parts.append(f"{symbol}, {formula}")
    
    print("; ".join(output_parts))

solve_polyhedra()