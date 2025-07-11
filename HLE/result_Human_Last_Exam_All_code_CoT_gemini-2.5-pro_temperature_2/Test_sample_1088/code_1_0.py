def calculate_primitives_for_carbon():
    """
    Calculates and explains the number of primitive Gaussians for a Carbon atom
    in the 6-311G** basis set.
    """
    atom_symbol = "C"
    atom_name = "Carbon"

    # The number of functions for each angular momentum (l)
    # s (l=0) has 1 function
    # p (l=1) has 3 functions (px, py, pz)
    # d (l=2) has 5 functions (using spherical harmonics, the modern standard)
    deg = {'s': 1, 'p': 3, 'd': 5}

    print(f"Calculating the number of primitive Gaussians in a 6-311G** basis set for {atom_name} ({atom_symbol}).\n")
    print("The total is a sum of primitives for core, valence, and polarization functions.\n")

    # 1. Core Primitives ('6-')
    # The 1s core orbital is a contraction of 6 's' type primitives.
    core_prims = 6 * deg['s']
    print(f"Step 1: Core '1s' orbital from the '6-' prefix.")
    print(f"   - Calculation: 6 primitives * {deg['s']} (for an s-orbital) = {core_prims} primitives.\n")

    # 2. Valence Primitives ('-311G')
    # For heavy atoms, s and p orbitals in the same valence shell are built
    # from the same set of primitive exponents (sp-shells).
    # Inner valence shell from '-3..'
    val_inner_prims = (3 * deg['s']) + (3 * deg['p'])
    print(f"Step 2: Valence '2sp' orbitals from the '-311G' part.")
    print(f"   - Inner shell ('-3..'): (3 primitives * {deg['s']} s-func) + (3 primitives * {deg['p']} p-func) = {val_inner_prims} primitives.")
    
    # Middle valence shell from '..1.'
    val_mid_prims = (1 * deg['s']) + (1 * deg['p'])
    print(f"   - Middle shell ('..1.'): (1 primitive * {deg['s']} s-func) + (1 primitive * {deg['p']} p-func) = {val_mid_prims} primitives.")

    # Outer valence shell from '...1'
    val_outer_prims = (1 * deg['s']) + (1 * deg['p'])
    print(f"   - Outer shell ('...1'): (1 primitive * {deg['s']} s-func) + (1 primitive * {deg['p']} p-func) = {val_outer_prims} primitives.\n")

    # 3. Polarization Primitives ('**')
    # The first '*' adds d-functions to heavy atoms. It's one uncontracted shell.
    pol_prims = 1 * deg['d']
    print(f"Step 3: Polarization functions from the '*' symbol.")
    print(f"   - One set of d-functions is added, using 1 primitive for each of the {deg['d']} d-functions.")
    print(f"   - Calculation: 1 primitive * {deg['d']} (d-functions) = {pol_prims} primitives.\n")
    print("(The second '*' adds p-functions to Hydrogen atoms, so it doesn't apply here).\n")


    # 4. Final Sum
    total_primitives = core_prims + val_inner_prims + val_mid_prims + val_outer_prims + pol_prims
    
    print("------------------------------------------------------------------")
    print("FINAL CALCULATION:")
    print("The total number of primitives is the sum of all components.")
    # Here we print each number in the final equation as requested
    print(f"Total = {core_prims} (core) + {val_inner_prims} (inner val) + {val_mid_prims} (mid val) + {val_outer_prims} (outer val) + {pol_prims} (pol) = {total_primitives}")
    print("------------------------------------------------------------------")


# Execute the function to print the result.
calculate_primitives_for_carbon()

# The final answer for Carbon is 31. For Hydrogen it would be 8.
# Since the question is general, the most common interpretation refers to a heavy atom like Carbon.
print("\n<<<31>>>")