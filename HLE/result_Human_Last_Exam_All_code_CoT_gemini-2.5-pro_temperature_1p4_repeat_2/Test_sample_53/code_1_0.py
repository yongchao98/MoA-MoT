def trace_pie_to_me():
    """
    Traces a hypothetical PIE root through sound changes to its
    potential Middle English form, printing each step of the derivation.
    """
    # Step 1: Proto-Indo-European (PIE)
    # The task asks for an o-grade causative. The o-grade of *kʷeys- is *kʷoys-.
    # The causative suffix is *-éye-. The 3rd person singular present indicative form is *kʷoyséyeti.
    pie_form = "*kʷoyséyeti"
    print("1. The Proto-Indo-European o-grade causative form is constructed.")
    print(f"   PIE form: {pie_form} (meaning 'he causes to heed/see', i.e., 'he shows')\n")

    # Step 2: Proto-Germanic (PGmc)
    # Sound changes from PIE to PGmc:
    # - Grimm's Law: *kʷ > *hʷ
    # - Verner's Law: *s > *z (because the PIE accent was on the suffix, not the root)
    # - Vowel Shift: *o > *a
    # - Suffix development: *-éyeti becomes the Class 1 weak verb stem *-azja- + endings.
    # The resulting infinitive is *hʷazjaną.
    pgmc_infinitive = "*hʷazjaną"
    print("2. The form undergoes changes into Proto-Germanic:")
    print("   - Grimm's Law (*kʷ > *hʷ) and Verner's Law (*s > *z) are applied.")
    print("   - Vowel and suffix changes result in a Class 1 weak verb.")
    print(f"   Proto-Germanic infinitive: {pgmc_infinitive}\n")

    # Step 3: Old English (OE)
    # Sound changes from PGmc to OE:
    # - Rhotacism: *z > *r in West Germanic. The stem becomes *hwarj-.
    # - I-umlaut: The *j fronts the stem vowel *a to *e. The stem becomes *hweri-.
    # - The 3rd person singular present ending for this class is -eþ.
    oe_form = "hwereþ"
    print("3. The form is inherited into Old English:")
    print("   - Rhotacism changes *z to r.")
    print("   - I-umlaut changes the stem vowel *a to e.")
    print("   - The 3rd person singular present ending '-eþ' is added.")
    print(f"   Old English form: {oe_form}\n")

    # Step 4: Middle English (ME)
    # Sound changes from OE to ME:
    # - Spelling change: 'hw' becomes 'wh'.
    # - Inflection: The ending '-eþ' is preserved and spelled '-eth'.
    me_form = "whereth"
    print("4. The form evolves into Middle English:")
    print("   - The spelling of initial 'hw' changes to 'wh'.")
    print("   - The ending '-eþ' becomes '-eth'.")
    print(f"   Middle English form: {me_form}\n")
    
    # Final "equation" showing the derivation path
    print("The complete derivation path is:")
    print(f"{pie_form} > {pgmc_infinitive} > {oe_form} > {me_form}")

if __name__ == '__main__':
    trace_pie_to_me()