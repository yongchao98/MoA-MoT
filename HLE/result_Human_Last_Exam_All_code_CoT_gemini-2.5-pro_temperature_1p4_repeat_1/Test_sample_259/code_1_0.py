import numpy as np
from collections import Counter

def find_coordination_polyhedra():
    """
    Analyzes the ReAl12 crystal structure to determine the coordination polyhedra for each unique atom site.
    The method involves generating a supercell, calculating interatomic distances, and identifying the first
    coordination shell by detecting a significant gap in distances.
    """

    # Step 1: Define crystal structure data from the problem description.
    # Cell parameters for the orthorhombic system
    cell_lengths = {
        'a': 7.609100,
        'b': 6.611700,
        'c': 9.023000
    }
    # The transformation matrix from fractional to Cartesian coordinates is diagonal for orthorhombic cells.
    cell_matrix = np.diag([cell_lengths['a'], cell_lengths['b'], cell_lengths['c']])

    # Atoms in the asymmetric unit
    asymmetric_unit_atoms = [
        {'label': 'Al_A', 'element': 'Al', 'coords': np.array([0.3182, 0.2158, 0.2500])},
        {'label': 'Al_B', 'element': 'Al', 'coords': np.array([0.0000, 0.3662, 0.1030])},
        {'label': 'Al_C', 'element': 'Al', 'coords': np.array([0.1743, 0.0000, 0.0000])},
        {'label': 'Re_A', 'element': 'Re', 'coords': np.array([0.0000, 0.0445, 0.2500])},
    ]

    # Symmetry operations for space group Cmcm (No. 63)
    symmetry_operations = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2', '-x+1/2,y+1/2,z',
        '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2', 'x+1/2,y+1/2,z',
        'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
    ]

    # Step 2: Generate all atom positions in the unit cell.
    unit_cell_atoms = []
    seen_coords_set = set()
    for atom_info in asymmetric_unit_atoms:
        for op_str in symmetry_operations:
            x, y, z = atom_info['coords']
            # Evaluate the symmetry operation string to get new coordinates.
            new_coords = np.array(eval(op_str))
            
            # Normalize coordinates to be within the primary unit cell [0, 1).
            new_coords = np.mod(new_coords, 1.0)
            
            # Use a rounded tuple of coordinates as a key to check for duplicates.
            coord_key = tuple(np.round(new_coords, 4))
            
            if coord_key not in seen_coords_set:
                seen_coords_set.add(coord_key)
                unit_cell_atoms.append({'element': atom_info['element'], 'coords': new_coords})

    # Step 3: Construct a 3x3x3 supercell to handle periodic boundaries.
    supercell_atoms = []
    for atom_info in unit_cell_atoms:
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    offset_vector = np.array([i, j, k])
                    supercell_atoms.append({
                        'element': atom_info['element'],
                        'coords': atom_info['coords'] + offset_vector
                    })

    # Step 4: For each unique atom, find its coordination polyhedron.
    results = []
    for central_atom_info in asymmetric_unit_atoms:
        central_coords = central_atom_info['coords']
        
        # Calculate distances to all atoms in the supercell.
        distances = []
        for neighbor_info in supercell_atoms:
            delta_frac = neighbor_info['coords'] - central_coords
            delta_cart = np.dot(delta_frac, cell_matrix)
            dist = np.linalg.norm(delta_cart)
            
            # Ignore the distance to the atom itself (which is close to zero).
            if dist > 1e-4:
                distances.append({'dist': dist, 'element': neighbor_info['element']})
        
        # Sort neighbors by distance.
        distances.sort(key=lambda item: item['dist'])
        
        # Step 5: Identify the first coordination shell using a distance gap threshold.
        if not distances:
            continue
            
        coordination_shell = [distances[0]]
        gap_threshold = 0.3  # Anstroms
        for i in range(1, len(distances)):
            if (distances[i]['dist'] - distances[i-1]['dist']) > gap_threshold:
                break
            coordination_shell.append(distances[i])
        
        # Step 6: Determine the chemical formula of the polyhedron.
        neighbor_elements = [item['element'] for item in coordination_shell]
        counts = Counter(neighbor_elements)
        
        # Sort elements by name (e.g., 'Al' then 'Re') for a consistent formula.
        formula = "".join(f"{el}{count}" for el, count in sorted(counts.items()))
        
        results.append(f"{central_atom_info['element']}, {formula}")
    
    # Step 7: Format and print the final, unique results.
    unique_results = sorted(list(set(results)))
    print("; ".join(unique_results))

find_coordination_polyhedra()