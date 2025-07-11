import pandas as pd

def solve_chemistry_problem():
    """
    Analyzes the stability of four Iridium complexes to determine which have shorter lifetimes.

    The stability and lifetime of these cyclometalated Ir(III) complexes are primarily
    determined by the strength of the Iridium-Carbon (Ir-C) bond.
    Electron-withdrawing fluorine substituents on the phenyl ring strengthen this bond,
    leading to greater stability and longer lifetimes.

    The stabilizing effect of fluorine depends on its position:
    - An ortho-fluorine (position 2') has a strong inductive withdrawing effect.
    - A para-fluorine (position 4') has a weaker inductive effect.

    A stability score is assigned based on this principle to rank the complexes.
    """
    # Define the complexes based on their cyclometalating (C^N) ligands
    # Each complex has two identical C^N ligands.
    complexes_data = {
        'Complex 1': {'ligand': '2-(2\',4\'-difluorophenyl)pyridine', 'fluorines': ['ortho', 'para']},
        'Complex 2': {'ligand': '2-phenylpyridine', 'fluorines': []},
        'Complex 3': {'ligand': '2-(2\'-fluorophenyl)pyridine', 'fluorines': ['ortho']},
        'Complex 4': {'ligand': '2-(4\'-fluorophenyl)pyridine', 'fluorines': ['para']}
    }

    # Define stability scores based on fluorine position
    stability_points = {
        'ortho': 2,
        'para': 1
    }

    # Calculate stability score for each complex
    results = []
    print("Step 1: Calculating stability scores based on fluorination.")
    print("Scoring rule: ortho-F = +2 points, para-F = +1 point per ligand.\n")

    for name, data in complexes_data.items():
        ligand_score = sum(stability_points.get(pos, 0) for pos in data['fluorines'])
        # Each complex has two such ligands
        total_score = 2 * ligand_score
        results.append({'Complex': name, 'Ligand': data['ligand'], 'Score': total_score})
        print(f"{name} ({data['ligand']}):")
        print(f"  - Score per ligand = {ligand_score}")
        print(f"  - Total Stability Score = 2 * {ligand_score} = {total_score}\n")

    # Create a pandas DataFrame for easy sorting and display
    df = pd.DataFrame(results)

    # Sort complexes by stability score in ascending order (least stable first)
    df_sorted = df.sort_values(by='Score', ascending=True)

    print("Step 2: Ranking complexes by stability (lower score = shorter lifetime).")
    print(df_sorted.to_string(index=False))
    print("\n")

    # The complexes with shorter lifetimes are the ones with the lowest scores.
    # We will identify the two least stable complexes.
    shorter_lifetime_complexes = df_sorted.head(2)
    
    print("Step 3: Identifying the complexes with the shortest lifetimes.")
    print("The two complexes with the lowest stability scores are expected to have shorter lifetimes.")
    
    complex_numbers = [c.split(' ')[1] for c in shorter_lifetime_complexes['Complex']]
    
    print(f"These are Complex {complex_numbers[0]} and Complex {complex_numbers[1]}.")

solve_chemistry_problem()
<<<I>>>