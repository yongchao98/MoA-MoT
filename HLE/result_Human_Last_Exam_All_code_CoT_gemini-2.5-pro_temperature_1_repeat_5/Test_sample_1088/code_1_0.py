import re

def calculate_primitives_for_formula(formula: str):
    """
    Calculates and explains the number of primitive Gaussians in a 6-311G**
    basis set for a given chemical formula.
    """
    print(f"Calculating primitive Gaussians for '{formula}' with the 6-311G** basis set.\n")

    # Data for atoms: number of core and valence atomic orbitals (AOs).
    # This data covers common main-group elements.
    # - Period 1: H
    # - Period 2: C, N, O (Core: 1s -> 1 AO; Valence: 2s, 2p -> 4 AOs)
    # - Period 3: S (Core: 1s, 2s, 2p -> 5 AOs; Valence: 3s, 3p -> 4 AOs)
    ATOMIC_DATA = {
        'H': {'name': 'Hydrogen', 'core_ao': 0, 'valence_ao': 1},
        'C': {'name': 'Carbon', 'core_ao': 1, 'valence_ao': 4},
        'N': {'name': 'Nitrogen', 'core_ao': 1, 'valence_ao': 4},
        'O': {'name': 'Oxygen', 'core_ao': 1, 'valence_ao': 4},
        'S': {'name': 'Sulfur', 'core_ao': 5, 'valence_ao': 4},
    }

    # This regex parses chemical formulas like "H2O" or "CH4"
    pattern = re.compile(r'([A-Z][a-z]?)(\d*)')
    matches = pattern.findall(formula)
    
    atom_counts = {}
    for symbol, count_str in matches:
        if symbol not in ATOMIC_DATA:
            print(f"Warning: Data for atom '{symbol}' is not available. Skipping this element.")
            continue
        count = int(count_str) if count_str else 1
        atom_counts[symbol] = atom_counts.get(symbol, 0) + count

    total_primitives = 0
    final_equation_parts = []

    print("--- Breakdown by Atom Type ---")

    # Process atoms alphabetically for consistent output
    for symbol, count in sorted(atom_counts.items()):
        data = ATOMIC_DATA[symbol]
        name = data['name']
        
        if symbol == 'H':
            # For Hydrogen: 5 from 311G valence + 3 from p-polarization
            valence_pgfs = 5
            polarization_pgfs = 3
            primitives_per_atom = valence_pgfs + polarization_pgfs
            
            print(f"Atom: {name} (H)")
            print(f"  Valence PGFs (from 311G): 3 + 1 + 1 = {valence_pgfs}")
            print(f"  Polarization PGFs (p-type from **): {polarization_pgfs}")
            print(f"  Total per H atom = {valence_pgfs} + {polarization_pgfs} = {primitives_per_atom}")
        else:
            # For heavy atoms
            core_ao = data['core_ao']
            valence_ao = data['valence_ao']
            
            # Core: 6 PGFs per core AO
            core_pgfs = core_ao * 6
            # Valence: 5 PGFs (3+1+1) per valence AO
            valence_pgfs = valence_ao * 5
            # Polarization: 6 d-type PGFs
            polarization_pgfs = 6
            primitives_per_atom = core_pgfs + valence_pgfs + polarization_pgfs

            print(f"Atom: {name} ({symbol})")
            print(f"  Core PGFs ({core_ao} AO(s) * 6): {core_pgfs}")
            print(f"  Valence PGFs ({valence_ao} AO(s) * 5): {valence_pgfs}")
            print(f"  Polarization PGFs (d-type from *): {polarization_pgfs}")
            print(f"  Total per {symbol} atom = {core_pgfs} + {valence_pgfs} + {polarization_pgfs} = {primitives_per_atom}")
        
        atom_total = primitives_per_atom * count
        total_primitives += atom_total
        final_equation_parts.append(f"({count} * {primitives_per_atom})") # (count * PGFs_per_atom)
        print(f"  Contribution for {count} atom(s): {count} * {primitives_per_atom} = {atom_total}\n")

    print("--- Final Calculation ---")
    final_equation = " + ".join(final_equation_parts)
    print(f"Total Primitive Gaussians = {final_equation} = {total_primitives}")


# --- Main Execution ---
# You can change this formula to calculate for other molecules like "H2O", "H2SO4", etc.
example_formula = "CH4"
calculate_primitives_for_formula(example_formula)