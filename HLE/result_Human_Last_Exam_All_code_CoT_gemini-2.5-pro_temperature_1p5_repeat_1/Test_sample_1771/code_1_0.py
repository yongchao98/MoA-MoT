def trace_etymology():
    """
    Traces the etymology of a PIE root through Germanic languages
    and prints the step-by-step derivation.
    """
    pie_form = "*seh₂gieti"
    pie_meaning = "'(s)he’s giving a sign'"

    # While the prompt gives the e-grade *seh₂gieti, the common Germanic verb 'to seek'
    # descends from the related o-grade causative PIE form, *soh₂géyeti 'to cause to track, to seek'.
    # We will trace this latter form, as it leads to the attested words.
    pie_causative = "*soh₂géyeti"
    pgmc_form = "*sōkiþi"

    # Stage 1: Proto-Indo-European to Proto-Germanic
    print(f"Tracing the development from PIE {pie_form} ({pie_meaning}).")
    print("The attested Germanic verb 'to seek' derives from the related PIE causative form.\n")
    print("1. From Proto-Indo-European to Proto-Germanic:")
    print(f"   PIE {pie_causative}  >  PGmc. {pgmc_form} '(s)he seeks'")
    print("   Sound Changes:")
    print("   - Vowel change: The PIE combination of a short vowel and a laryngeal (*o + *h₂) becomes a long vowel (*ō).")
    print("   - Grimm's Law: The voiced stop *g shifts to a voiceless stop *k.")
    print("   - Verb ending: The PIE 3rd singular ending *-eti (part of the causative *-éyeti suffix) becomes *-iþi in Proto-Germanic for this verb class.\n")

    # Stage 2: Proto-Germanic to Daughter Languages
    on_form = "sœkir"
    oe_form = "sēċeþ"
    ohg_form = "suohhit"

    print(f"2. From Proto-Germanic {pgmc_form} to its daughter languages:\n")

    # Old Norse
    print("   A. Old Norse:")
    print(f"      PGmc. {pgmc_form}  >  ON {on_form}")
    print("      - I-Umlaut: The back vowel *ō becomes the front vowel œ when an *i or *j follows in the next syllable.")
    print("      - Ending: The 3rd person singular present ending *-iþi becomes -r.\n")

    # Old English
    print("   B. Old English:")
    print(f"      PGmc. {pgmc_form}  >  OE {oe_form}")
    print("      - I-Umlaut: The long vowel *ō is fronted to *œ̄, which then becomes ē in Old English.")
    print("      - Palatalization: The consonant *k becomes the affricate ċ (pronounced /tʃ/) before a front vowel (like the ē that resulted from umlaut).")
    print("      - Ending: The ending *-iþi becomes -eþ.\n")

    # Old High German
    print("   C. Old High German:")
    print(f"      PGmc. {pgmc_form}  >  OHG {suohhit}")
    print("      - High German Diphthongization: The long vowel *ō becomes the diphthong uo.")
    print("      - Second Sound Shift: The voiceless stop *k becomes a double fricative hh /xx/ after a vowel.")
    print("      - Ending: The ending *-iþi becomes -it.\n")

    # Final summary
    print("="*40)
    print("Summary of the derived forms:")
    print(f"Proto-Germanic: {pgmc_form}")
    print(f"Old Norse: {on_form}")
    print(f"Old English: {oe_form}")
    print(f"Old High German: {ohg_form}")
    print("="*40)


if __name__ == '__main__':
    trace_etymology()
<<<Proto-Germanic: *sōkiþi, Old Norse: sœkir, Old English: sēċeþ, Old High German: suohhit>>>