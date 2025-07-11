def analyze_complex_stability():
    """
    Analyzes the stability of four Iridium complexes based on their structure.
    The lifetime of these complexes in LECs is related to their stability.
    A key factor for instability is steric hindrance caused by substituents
    at the ortho-position to the Ir-C bond on the phenyl ring.
    """

    # Data represents the presence (True) or absence (False) of a fluorine
    # substituent at the critical ortho-position for each complex.
    complexes_info = {
        1: {'has_ortho_F': False, 'fluorination': '4,5-difluoro'},
        2: {'has_ortho_F': True, 'fluorination': '2,4-difluoro'},
        3: {'has_ortho_F': True, 'fluorination': '2-fluoro'},
        4: {'has_ortho_F': True, 'fluorination': '2,3,4-trifluoro'}
    }

    shorter_lifetime_complexes = []
    print("Analysis of Complex Stability and Expected Lifetime:")
    print("-" * 50)
    print("The stability of these Ir(III) complexes is a key predictor of their lifetime in an LEC.")
    print("A major cause of degradation is ligand dissociation, which can be induced by steric strain.\n")
    print("The critical structural feature is a substituent at the ortho-position of the phenyl ring (relative to the Ir-C bond). This causes steric clash with the bipyridine ligand, weakening the complex.\n")

    for i in range(1, 5):
        info = complexes_info[i]
        ortho_status = "present" if info['has_ortho_F'] else "absent"
        stability_pred = "less stable (shorter lifetime)" if info['has_ortho_F'] else "more stable (longer lifetime)"
        print(f"Complex {i} ({info['fluorination']}):")
        print(f"  - Ortho-fluorine substituent: {ortho_status}")
        print(f"  - Expected stability: {stability_pred}")
        if info['has_ortho_F']:
            shorter_lifetime_complexes.append(i)
    
    print("-" * 50)
    print("Conclusion: The complexes with an ortho-fluorine substituent are expected to have shorter lifetimes due to increased steric strain.")
    print("The complexes predicted to have shorter lifetimes are:")
    
    # Format the output as requested, e.g., [2, 3, 4]
    # The prompt asks to "output each number in the final equation"
    final_list_str = f"[{', '.join(map(str, shorter_lifetime_complexes))}]"
    print(final_list_str)

analyze_complex_stability()