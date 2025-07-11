def calculate_yields():
    """
    Analyzes the reaction, calculates molar masses, identifies limiting reagents,
    and computes the theoretical yields for products A, B, and C.
    """
    # Step 1: Define constants and molecular formulas
    atomic_weights = {'C': 12.01, 'H': 1.008, 'N': 14.01, 'O': 16.00}
    
    # Reactant properties
    sm_mass = 1.0  # g
    mp_volume = 4.0  # mL
    mp_density = 1.08  # g/mL
    ac2o_volume = 10.0  # mL
    ac2o_density = 1.082 # g/mL

    # Formulas from analysis
    formulas = {
        'SM': {'C': 9, 'H': 14, 'N': 2, 'O': 2},
        'MP': {'C': 4, 'H': 4, 'O': 2},
        'Ac2O': {'C': 4, 'H': 6, 'O': 3},
        'A': {'C': 14, 'H': 20, 'N': 2, 'O': 3},
        'B': {'C': 12, 'H': 14, 'N': 2, 'O': 3},
        'C': {'C': 11, 'H': 16, 'N': 2, 'O': 3},
    }

    # Step 2: Define a function to calculate Molar Mass (MW)
    def get_mw(molecule_name):
        formula = formulas[molecule_name]
        mw = 0
        calculation_str = f"MW({molecule_name}) = "
        terms = []
        for atom, count in formula.items():
            weight = atomic_weights[atom]
            mw += count * weight
            terms.append(f"{count} * {weight}")
        calculation_str += " + ".join(terms)
        calculation_str += f" = {mw:.2f} g/mol"
        print(calculation_str)
        return mw

    print("--- Molar Mass Calculations ---")
    mw_sm = get_mw('SM')
    mw_mp = get_mw('MP')
    mw_ac2o = get_mw('Ac2O')
    mw_a = get_mw('A')
    mw_b = get_mw('B')
    mw_c = get_mw('C')
    print("-" * 30)

    # Step 3: Calculate initial moles of reactants
    moles_sm = sm_mass / mw_sm
    
    mass_mp = mp_volume * mp_density
    moles_mp = mass_mp / mw_mp
    
    mass_ac2o = ac2o_volume * ac2o_density
    moles_ac2o = mass_ac2o / mw_ac2o
    
    print("\n--- Initial Moles of Reactants ---")
    print(f"Starting Material (SM): {moles_sm:.4f} mol")
    print(f"Methyl Propiolate (MP): {moles_mp:.4f} mol")
    print(f"Acetic Anhydride (Ac2O): {moles_ac2o:.4f} mol")
    print("-" * 30)

    # Step 4: Calculate Theoretical Yield for each product
    print("\n--- Theoretical Yield Calculations ---")
    
    # For Product A (and C): Reaction is 1 SM + 1 Ac2O + 1 MP -> 1 A ...
    # Limiting reagent is SM since 0.0055 < 0.1060 (Ac2O) and 0.0055 < 0.0514 (MP)
    moles_a_theoretical = moles_sm
    yield_a = moles_a_theoretical * mw_a
    print(f"Product A is formed from SM, Ac2O, and MP. Limiting reagent is SM.")
    print(f"Theoretical Yield of A: {moles_a_theoretical:.4f} mol * {mw_a:.2f} g/mol = {yield_a:.2f} g")

    # For Product B: Reaction is 1 SM + 1 MP -> 1 B ...
    # Limiting reagent is SM since 0.0055 < 0.0514 (MP)
    moles_b_theoretical = moles_sm
    yield_b = moles_b_theoretical * mw_b
    print(f"Product B is formed from SM and MP. Limiting reagent is SM.")
    print(f"Theoretical Yield of B: {moles_b_theoretical:.4f} mol * {mw_b:.2f} g/mol = {yield_b:.2f} g")

    # For Product C: Reaction is 1 SM + 1 Ac2O -> 1 C ...
    # Limiting reagent is SM since 0.0055 < 0.1060 (Ac2O)
    moles_c_theoretical = moles_sm
    yield_c = moles_c_theoretical * mw_c
    print(f"Product C is formed from SM and Ac2O. Limiting reagent is SM.")
    print(f"Theoretical Yield of C: {moles_c_theoretical:.4f} mol * {mw_c:.2f} g/mol = {yield_c:.2f} g")
    print("-" * 30)

calculate_yields()