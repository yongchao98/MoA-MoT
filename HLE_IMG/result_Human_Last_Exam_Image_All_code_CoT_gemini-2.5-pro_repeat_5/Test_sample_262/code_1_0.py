def find_short_lifetime_complexes():
    """
    Identifies iridium complexes expected to have shorter lifetimes based on ligand substitution.

    The stability of these Ir(III) complexes in LECs is largely dictated by the strength of the
    cyclometalating Ir-C bond. Electron-withdrawing groups, such as fluorine, on the phenyl
    ring weaken this bond, leading to faster degradation and shorter device lifetimes.

    The position of the fluorine is critical. Fluorination at the 2'-position (ortho to the
    Ir-C bond) is known to be particularly destabilizing.

    This program will identify which of the four complexes contain this 2'-fluoro feature.
    """
    complexes = {
        1: {'name': 'Complex 1', 'substituents': ['2-F', '4-F']},
        2: {'name': 'Complex 2', 'substituents': []},
        3: {'name': 'Complex 3', 'substituents': ['2-F']},
        4: {'name': 'Complex 4', 'substituents': ['4-F']}
    }

    short_lifetime_complex_ids = []
    print("Analysis of Complexes for Lifetime Stability:")
    print("-" * 45)
    print("Key Destabilizing Feature: Fluorine at the 2'-position (ortho-position).\n")

    for i, data in complexes.items():
        if '2-F' in data['substituents']:
            short_lifetime_complex_ids.append(i)
            print(f"- {data['name']} [ID: {i}]: Contains a 2'-fluoro substituent. Expected to have a SHORT lifetime.")
        else:
            if not data['substituents']:
                print(f"- {data['name']} [ID: {i}]: No fluoro substituents. Expected to have the LONGEST lifetime.")
            else:
                 print(f"- {data['name']} [ID: {i}]: Contains a 4'-fluoro substituent. Lifetime is shorter than unsubstituted but longer than 2'-fluorinated complexes.")


    print("\nConclusion:")
    print("Complexes with the 2'-fluoro substituent are the most unstable.")
    print(f"The complexes expected to show shorter lifetimes are: {sorted(short_lifetime_complex_ids)}")

find_short_lifetime_complexes()