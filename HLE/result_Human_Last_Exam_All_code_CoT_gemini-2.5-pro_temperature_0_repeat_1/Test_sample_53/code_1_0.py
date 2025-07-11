def solve_linguistic_derivation():
    """
    Traces the hypothetical derivation of a PIE root to a Middle English verb form
    and prints the step-by-step evolution.
    """
    # Define the forms at each stage of the language evolution
    pie_root = "*kʷeys-"
    pgmc_form = "*hʷaizijaną"
    pwgmc_form = "*hʷaizan"
    oe_form = "hwāreþ"
    me_form = "whoreth"

    # Print the derivation explanation
    print("This script traces the hypothetical evolution of a PIE root into a Middle English verb.")
    print("The derivation follows standard historical sound changes for an o-grade causative verb.")
    print("\n--- Derivation Path ---\n")

    print(f"1. Proto-Indo-European Root: {pie_root} ('to see, to heed')")
    print(f"   - Changes to Proto-Germanic:")
    print(f"     a. O-grade (*e > *o), Grimm's Law (*kʷ > *hʷ), Verner's Law (*s > *z), Vowel Shift (*oy > *ai)")
    print(f"   ---> Proto-Germanic: {pgmc_form}\n")

    print(f"2. Proto-Germanic -> Proto-West Germanic:")
    print(f"   - Changes:")
    print(f"     a. Loss of *j after a long syllable.")
    print(f"   ---> Proto-West Germanic: {pwgmc_form}\n")

    print(f"3. Proto-West Germanic -> Old English:")
    print(f"   - Changes:")
    print(f"     a. *hʷ > hw, *ai > ā, *z > r (rhotacism)")
    print(f"     b. 3rd person singular present ending '-eþ' is added.")
    print(f"   ---> Old English: {oe_form}\n")

    print(f"4. Old English -> Middle English:")
    print(f"   - Changes:")
    print(f"     a. Spelling hw > wh, Vowel Shift ā > ō (spelled 'o'), ending -eþ > -eth")
    print(f"   ---> Middle English: {me_form}\n")

    # Print the final summary equation as requested
    print("--- Final Equation ---")
    final_equation = f"{pie_root} > {pgmc_form} > {pwgmc_form} > {oe_form} > {me_form}"
    print(final_equation)

solve_linguistic_derivation()
<<<whoreth>>>