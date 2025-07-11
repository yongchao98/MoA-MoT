def derive_middle_english_form():
    """
    This function traces a hypothetical PIE root through its linguistic evolution
    to a predicted Middle English form, printing each step of the derivation.
    """

    print("This script will derive a hypothetical Middle English word from a Proto-Indo-European (PIE) root.")
    print("------------------------------------------------------------------------------------------------\n")

    # --- Stage 1: Proto-Indo-European (PIE) ---
    print("Step 1: Construct the Proto-Indo-European Form")
    pie_root = "*kʷeys-"
    pie_o_grade_stem = "*kʷoys-"
    pie_3sg_pres_causative = "*kʷoyséyeti"
    print(f"The starting root is {pie_root} (to see, to heed).")
    print(f"The o-grade causative stem is formed by changing the root vowel to 'o', resulting in '{pie_o_grade_stem}'.")
    print("The meaning becomes 'to cause to see', i.e., 'to show'.")
    print(f"The 3rd person singular present indicative form is: {pie_3sg_pres_causative}\n")

    # --- Stage 2: Proto-Germanic (PGmc) ---
    print("Step 2: Apply Sound Changes to form Proto-Germanic")
    pgmc_3sg_pres = "*hʷaiziþi"
    print(f" - Grimm's Law: PIE *kʷ becomes PGmc *hʷ.")
    print(f" - Verner's Law: Because the PIE accent was on the suffix, the preceding 's' becomes voiced, giving PGmc *z.")
    print(f" - Vowel Shift: The PIE diphthong *oi becomes PGmc *ai.")
    print(f" - Suffix Change: The verb enters Class 1 of weak verbs, and the 3sg. ending becomes *-iþi.")
    print(f"The resulting Proto-Germanic form is: {pgmc_3sg_pres}\n")

    # --- Stage 3: Old English (OE) ---
    print("Step 3: Apply Sound Changes to form Old English")
    oe_3sg_pres = "hwǣreþ"
    print(f" - Consonant Shift: PGmc *hʷ becomes 'hw' in Old English.")
    print(f" - i-Umlaut: PGmc *ai becomes OE 'ā', which is then fronted to 'ǣ' by the 'i' in the suffix.")
    print(f" - Rhotacism: Intervocalic *z becomes 'r' in Old English.")
    print(f" - Suffix Change: The PGmc ending *-iþi becomes OE '-eþ'.")
    print(f"The resulting Old English form is: {oe_3sg_pres}\n")

    # --- Stage 4: Middle English (ME) ---
    print("Step 4: Apply Sound Changes to form Middle English")
    me_3sg_pres = "whereth"
    print(f" - Spelling Change: The cluster 'hw' is commonly spelled 'wh' in Middle English.")
    print(f" - Vowel Shift: The long vowel OE 'ǣ' [æː] is raised to 'ē' [ɛː], which is typically spelled with 'e'.")
    print(f" - Suffix Simplification: The OE ending '-eþ' is preserved as '-eth'.")
    print(f"The resulting Middle English form is: {me_3sg_pres}\n")

    # --- Final Answer ---
    print("------------------------------------------------------------------------------------------------")
    print(f"The final predicted third person singular present form in Middle English is '{me_3sg_pres}'.")


# Execute the derivation function
derive_middle_english_form()
<<<whereth>>>