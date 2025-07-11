import sys

def solve_etymology():
    """
    This function traces a hypothetical PIE root to a Middle English verb form
    and prints the derivation process.
    """
    
    # Disable writing bytecode to avoid __pycache__ folders
    sys.dont_write_bytecode = True

    print("Deriving the Middle English third-person singular present form of a hypothetical verb from PIE *kʷeys-.\n")

    print("Step 1: Proto-Indo-European to Proto-Germanic")
    print("---------------------------------------------")
    print("1. Start with the PIE root: *kʷeys- (to see, to heed)")
    print("2. Form the o-grade causative stem: *kʷoys-eye- ('to make see', 'to show')")
    print("3. Apply Proto-Germanic sound changes:")
    print("   - Grimm's Law: PIE *kʷ > PGmc *hw")
    print("   - Vowel Shift: PIE *o > PGmc *a")
    print("   - Suffix Change: The causative suffix *-eye- forms a Class I weak verb, resulting in the infinitive suffix *-janą.")
    print("   => The Proto-Germanic infinitive is *hwaisijaną.")
    print("\n")

    print("Step 2: Proto-Germanic to Old English")
    print("---------------------------------------")
    print("1. The third-person singular form is derived from the PGmc stem *hwaisija-.")
    print("2. The 'j' in the stem causes i-umlaut, changing the root vowel 'ai' to 'ǣ' in Old English.")
    print("3. The initial PGmc consonant cluster *hw- is preserved as 'hw'.")
    print("4. The third-person singular present ending for Class I weak verbs becomes '-eþ'.")
    print("   => The Old English form is hwǣseþ.")
    print("\n")
    
    print("Step 3: Old English to Middle English")
    print("--------------------------------------")
    print("1. Apply sound and spelling changes from Old to Middle English:")
    print("   - Initial 'hw' becomes 'wh'.")
    print("   - The long vowel 'ǣ' /æː/ raises and becomes the long open 'ē' /ɛː/, spelled 'e'.")
    print("   - The ending '-eþ' is preserved and written as '-eth'.")
    print("   => The Middle English form is 'wheseth'.")
    print("\n")

    print("Final Derivation Summary:")
    print("PIE *kʷoys-éye-(ti) > PGmc *hwaisiþ(i) > OE hwǣseþ > ME wheseth\n")
    
    print("Final Middle English Form (meaning 'he shows'):")
    part1 = "wh"
    part2 = "e"
    part3 = "s"
    part4 = "eth"
    result = part1 + part2 + part3 + part4
    print(f"{part1} + {part2} + {part3} + {part4} = {result}")

solve_etymology()
<<<wheseth>>>