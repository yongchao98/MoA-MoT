import textwrap

def analyze_lec_stability():
    """
    Analyzes the stability of three Ir(III) complexes for use in LECs.
    """

    analysis = {
        "Complex 1": {
            "Structure": "[Ir(ppy)2(bpy)]+ (ppy = 2-phenylpyridine, bpy = 2,2'-bipyridine)",
            "Analysis": "This is a standard, widely-used cationic Ir(III) complex. It serves as a baseline for comparison. It lacks specific modifications for enhanced stability.",
            "Stability Factor": "Baseline"
        },
        "Complex 2": {
            "Structure": "[Ir(ppy)2(L)]+, where L is a large, sterically bulky N^N ligand.",
            "Analysis": "This complex is modified on the ancillary N^N ligand. It features large phenyl and tolyl groups. These bulky groups provide significant steric hindrance, which can prevent the complexes from packing too closely together. This reduces concentration quenching and aggregation-induced degradation, leading to improved device stability compared to Complex 1.",
            "Stability Factor": "Increased due to steric hindrance."
        },
        "Complex 3": {
            "Structure": "[Ir(dfppy)2(dtbbpy)]+ (dfppy = 2-(2,4-difluorophenyl)pyridine, dtbbpy = 4,4'-di-tert-butyl-2,2'-bipyridine)",
            "Analysis": textwrap.fill(
                "This complex incorporates two key strategies for enhancing stability. "
                "First, the cyclometalating (ppy) ligands are fluorinated. The strong electron-withdrawing effect of fluorine lowers the energy of the Highest Occupied Molecular Orbital (HOMO), making the complex more resistant to oxidative degradation. "
                "Second, the bipyridine ligand has bulky tert-butyl groups (dtbbpy). Similar to Complex 2, these groups provide steric protection against aggregation. The combination of electronic stabilization and steric hindrance makes this complex design very robust.",
                width=80
            ),
            "Stability Factor": "Greatly increased due to both steric hindrance and electronic stabilization (fluorination)."
        }
    }

    print("--- Analysis of LEC Emitter Stability ---")
    for complex_name, data in analysis.items():
        print(f"\nAnalyzing {complex_name}:")
        print(f"  - Structure Summary: {data['Structure']}")
        print(f"  - Stability Analysis: {data['Analysis']}")
        print(f"  - Resulting Stability: {data['Stability Factor']}")

    print("\n--- Conclusion ---")
    conclusion_text = textwrap.fill(
        "Based on the analysis, both Complex 2 and Complex 3 incorporate well-known molecular design strategies to improve stability over the baseline Complex 1. "
        "Complex 2 uses steric bulk, while Complex 3 uses both steric bulk and electronic stabilization from fluorination. "
        "Therefore, LECs based on both Complex 2 and Complex 3 are expected to be more stable than those based on Complex 1. "
        "The most appropriate answer choice that reflects this is F.",
        width=80
    )
    print(conclusion_text)

    # Final Answer Choice
    final_answer = "F"
    print(f"\nThe best answer choice is: {final_answer}")

# Run the analysis
analyze_lec_stability()
<<<F>>>