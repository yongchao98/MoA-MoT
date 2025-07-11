def count_nmr_peaks():
    """
    Calculates the number of expected 1H NMR peaks for the given molecule based on its structure and symmetry.
    
    The molecule is: 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene.
    
    The logic is as follows:
    1. The molecule has C3 symmetry, so we only need to analyze one of the three identical 'arms'.
    2. We count the number of chemically non-equivalent proton environments in one arm and its associated methyl group on the central ring.
    """
    
    # A list of tuples, where each tuple contains a description of the proton environment
    # and the number of distinct NMR signals it produces.
    signal_sources = [
        ("Central ring methyl groups (-CH3)", 1),
        ("Linker methylene protons (-CH2-), diastereotopic", 2),
        ("Pyrazole ring proton (C3-H)", 1),
        ("Camphor skeleton bridgehead proton (C4-H)", 1),
        ("Camphor skeleton bridgehead methyl group (C7-CH3)", 1),
        ("Camphor skeleton gem-dimethyl protons (C8-Me2), diastereotopic", 2),
        ("Camphor skeleton methylene protons (C5-H2), diastereotopic", 2),
        ("Camphor skeleton methylene protons (C6-H2), diastereotopic", 2)
    ]
    
    print("Analysis of unique proton environments:")
    total_signals = 0
    equation_parts = []
    
    for description, count in signal_sources:
        print(f"- {description}: {count} signal(s)")
        total_signals += count
        equation_parts.append(str(count))
        
    print("\nTotal number of expected peaks is the sum of signals from each unique environment.")
    equation = " + ".join(equation_parts)
    print(f"Calculation: {equation} = {total_signals}")
    
    return total_signals

if __name__ == "__main__":
    count_nmr_peaks()