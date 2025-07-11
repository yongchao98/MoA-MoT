def format_formula(atoms_dict):
    """Formats a dictionary of atoms into a standard molecular formula string."""
    # Order: C, H, then alphabetical for the rest.
    order = ['C', 'H', 'F', 'N', 'O']
    formula = ""
    for element in order:
        if atoms_dict.get(element, 0) > 0:
            formula += element
            if atoms_dict[element] > 1:
                formula += str(atoms_dict[element])
    return formula

def solve_reaction_sequence():
    """Calculates the molecular formula of the product through a three-step reaction."""
    
    print("Step 1: Determine the molecular formula of the starting material.")
    # The starting material is a 2-azabicyclo[2.2.1]hept-5-en-3-one derivative.
    # Core (C6H5NO): 6 C, 5 H, 1 N, 1 O
    # CF3 substituent: 1 C, 3 F
    # PMB (p-methoxybenzyl) substituent (C8H9O): 8 C, 9 H, 1 O
    c_start = 6 + 1 + 8
    h_start = 5 + 9
    f_start = 3
    n_start = 1
    o_start = 1 + 1
    
    start_atoms = {'C': c_start, 'H': h_start, 'F': f_start, 'N': n_start, 'O': o_start}
    print(f"Starting Material Formula: {format_formula(start_atoms)}")
    print(f"Initial counts: C={c_start}, H={h_start}, F={f_start}, N={n_start}, O={o_start}\n")
    
    print("Step 2: First reaction (PMB deprotection).")
    print("This reaction removes the PMB group (C8H9O) and adds one Hydrogen atom.")
    c_inter1 = c_start - 8
    h_inter1 = h_start - 9 + 1
    o_inter1 = o_start - 1
    inter1_atoms = {'C': c_inter1, 'H': h_inter1, 'F': f_start, 'N': n_start, 'O': o_inter1}
    print(f"Intermediate 1 Formula: {format_formula(inter1_atoms)}")
    print(f"Atom counts: C={c_inter1}, H={h_inter1}, F={f_start}, N={n_start}, O={o_inter1}\n")

    print("Step 3: Second reaction (Hydrogenation).")
    print("This adds one molecule of H2 across the C=C bond.")
    h_inter2 = h_inter1 + 2
    inter2_atoms = {'C': c_inter1, 'H': h_inter2, 'F': f_start, 'N': n_start, 'O': o_inter1}
    print(f"Intermediate 2 Formula: {format_formula(inter2_atoms)}")
    print(f"Atom counts: C={c_inter1}, H={h_inter2}, F={f_start}, N={n_start}, O={o_inter1}\n")

    print("Step 4: Third reaction (Lactam hydrolysis).")
    print("This adds one molecule of H2O to open the ring.")
    h_final = h_inter2 + 2
    o_final = o_inter1 + 1
    final_atoms = {'C': c_inter1, 'H': h_final, 'F': f_start, 'N': n_start, 'O': o_final}
    print(f"Final Product Formula: {format_formula(final_atoms)}\n")

    print("-" * 40)
    print("Summary of atom count changes to find the final formula:")
    print(f"Carbon (C)   = {c_start} - 8 = {c_inter1}")
    print(f"Hydrogen (H) = {h_start} - 9 + 1 + 2 + 2 = {h_final}")
    print(f"Fluorine (F) = {f_start}")
    print(f"Nitrogen (N) = {n_start}")
    print(f"Oxygen (O)   = {o_start} - 1 + 1 = {o_final}")
    print("-" * 40)

    print("The final molecular formula equation is:")
    print(f"C{final_atoms['C']}H{final_atoms['H']}F{final_atoms['F']}N{final_atoms['N']}O{final_atoms['O']}")

# Execute the calculation
solve_reaction_sequence()