import numpy as np
from collections import Counter

def find_coordination_polyhedra():
    """
    This script determines the coordination polyhedra for a crystal structure
    based on its CIF data. It generates all atoms in the unit cell,
    calculates distances to find nearest neighbors, and reports the
    resulting polyhedra formulas.
    """
    # Step 1: Parse Crystal Data from the problem description
    cell_lengths = np.array([7.609100, 6.611700, 9.023000])
    
    # The unique atoms in the asymmetric unit
    unique_atoms = [
        {"label": "Al_A", "symbol": "Al", "coords": np.array([0.3182, 0.2158, 0.2500])},
        {"label": "Al_B", "symbol": "Al", "coords": np.array([0.0000, 0.3662, 0.1030])},
        {"label": "Al_C", "symbol": "Al", "coords": np.array([0.1743, 0.0000, 0.0000])},
        {"label": "Re_A", "symbol": "Re", "coords": np.array([0.0000, 0.0445, 0.2500])}
    ]

    # Symmetry operations for space group Cmcm (No. 63)
    sym_ops_str = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2',
        '-x+1/2,y+1/2,z', '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2',
        'x+1/2,y+1/2,z', 'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
    ]

    # Step 2: Generate the Full Unit Cell
    cell_atoms = []
    seen_coords = set()
    for atom in unique_atoms:
        x, y, z = atom["coords"]
        for op in sym_ops_str:
            # Safely evaluate the symmetry operation string
            new_coords = np.array(eval(op, {"x": x, "y": y, "z": z}))
            # Normalize coordinates to be within the unit cell [0, 1)
            new_coords = new_coords - np.floor(new_coords)
            
            # Use a formatted string to check for duplicate positions with a tolerance
            coord_key = tuple(np.round(new_coords, 4))
            if coord_key not in seen_coords:
                cell_atoms.append({"symbol": atom["symbol"], "coords": new_coords})
                seen_coords.add(coord_key)

    # Step 3: Construct a Supercell for neighbor search
    supercell_atoms = []
    for atom in cell_atoms:
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    supercell_atoms.append({
                        "symbol": atom["symbol"],
                        "coords": atom["coords"] + np.array([i, j, k])
                    })

    # Step 4-7: Determine coordination, format, and collect results
    final_results = []
    for center_atom in unique_atoms:
        center_coords = center_atom["coords"]
        
        # Calculate distances to all atoms in the supercell
        distances = []
        for neighbor in supercell_atoms:
            delta_frac = neighbor["coords"] - center_coords
            
            # Skip the central atom itself
            if np.allclose(delta_frac, [0, 0, 0]):
                continue
            
            # Convert fractional distance vector to Cartesian for orthorhombic cell
            delta_cart = delta_frac * cell_lengths
            dist = np.linalg.norm(delta_cart)
            distances.append({"dist": dist, "symbol": neighbor["symbol"]})
        
        # Sort neighbors by distance
        distances.sort(key=lambda d: d['dist'])
        
        # Identify the first coordination shell by finding a large gap in distances
        if not distances:
            continue
            
        first_shell = [distances[0]]
        last_dist = distances[0]['dist']
        # Use a 20% distance increase as the cutoff for the first shell
        for i in range(1, len(distances)):
            current_dist = distances[i]['dist']
            if current_dist > last_dist * 1.20:
                break
            first_shell.append(distances[i])
            # Update last_dist only with non-duplicate distances to handle degeneracies
            if not np.isclose(current_dist, last_dist):
                last_dist = current_dist

        # Count elements in the shell to create the formula
        shell_symbols = [d['symbol'] for d in first_shell]
        counts = Counter(shell_symbols)
        
        # Build the formula string, sorted alphabetically by element
        formula = "".join(
            f"{symbol}{counts[symbol] if counts[symbol] > 1 else ''}"
            for symbol in sorted(counts.keys())
        )
        
        final_results.append(f"{center_atom['label']}, {formula}")

    print("; ".join(final_results) + ";")

find_coordination_polyhedra()