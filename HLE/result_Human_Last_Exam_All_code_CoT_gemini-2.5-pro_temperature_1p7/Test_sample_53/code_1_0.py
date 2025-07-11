def solve_linguistic_derivation():
    """
    This function traces a hypothetical PIE root through its sound changes
    into Middle English, following the rules of historical linguistics.
    It prints each stage of the derivation.
    """

    # Stage 1: Proto-Indo-European (PIE)
    # The root is *kʷeys- ("to see, to heed").
    # To form an o-grade causative, we change the root vowel 'e' to 'o'
    # and add the causative suffix *-éye-.
    pie_causative = "*kʷoys-éye-"
    print(f"1. Proto-Indo-European (Causative): {pie_causative}")

    # Stage 2: Proto-Germanic (PGmc)
    # This involves several key sound changes from PIE:
    # - Grimm's Law: The voiceless stop *kʷ becomes the voiceless fricative *hʷ.
    # - Verner's Law: Because the PIE accent was on the suffix (*-éye-), the
    #   preceding voiceless fricative *s becomes voiced, resulting in *z.
    # - Vowel Change: The PIE diphthong *oy becomes PGmc *ai.
    # - Suffix Change: The PIE suffix *-éye- becomes the PGmc Class 1 weak
    #   verb suffix *-janą.
    pgmc_verb = "*hʷaizijaną"
    print(f"   Sound Laws: *kʷ > *hʷ; *s > *z; *oy > *ai; *-éye- > *-janą")
    print(f"2. Proto-Germanic: {pgmc_verb}")

    # Stage 3: Old English (OE)
    # From Proto-Germanic, the word enters Old English via Proto-West Germanic.
    # - Rhotacism: In West Germanic, intervocalic *z became *r.
    #   *hʷaizijaną > *hʷairijan
    # - I-umlaut: The *j in the suffix causes the back vowel *ai to front. The vowel
    #   *ai first becomes *ā before *r, which then fronts to *ǣ (/æː/) due to umlaut.
    #   The *j is subsequently lost.
    # - Other changes: PGmc *hʷ > OE hw.
    # The resulting OE infinitive is *hǣran ("to cause to see/heed, to show").
    oe_infinitive = "hǣran"
    print(f"   Sound Laws: *z > *r; i-umlaut (*ā > *ǣ); *hʷ > hw")
    print(f"3. Old English (Infinitive): {oe_infinitive}")

    # For the 3rd person singular present tense, Class 1 weak verbs add the ending -eþ.
    oe_3rd_person = "hǣreþ"
    print(f"4. Old English (3rd Person Singular Present): {oe_3rd_person}")

    # Stage 4: Middle English (ME)
    # The final step is from Old to Middle English.
    # - Vowel Shift: The OE long vowel *ǣ (/æː/) is raised and becomes the
    #   Middle English open *ē (/ɛː/), typically spelled 'e'.
    # - Spelling Change: The OE consonant cluster 'hw' becomes 'wh' in spelling.
    # - Inflectional Leveling: The OE ending -eþ becomes the standard ME ending -eth.
    me_3rd_person = "whereth"
    print(f"   Sound Laws: OE ǣ > ME ē; OE hw > ME wh; OE -eþ > ME -eth")
    print(f"5. Middle English (3rd Person Singular Present): {me_3rd_person}")

    # The final derived word form means "he/she shows".
    final_answer = me_3rd_person
    print("\n---")
    print(f"The final predicted Middle English form is: {final_answer}")
    return final_answer

solve_linguistic_derivation()