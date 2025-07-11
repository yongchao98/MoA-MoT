def derive_middle_english_form():
    """
    This function derives a hypothetical Middle English verb form
    by applying a series of historical sound changes.
    """
    print("Deriving the Middle English reflex of PIE *kʷeys- ('to see').")
    print("Target: 3rd person singular, present tense, o-grade causative ('he shows').")
    print("-" * 70)

    # Step 1: Proto-Indo-European (PIE)
    # O-grade of *kʷeys- is *kʷoys-. Causative is *kʷoys-éye-. 3sg pres is *kʷoys-éye-ti.
    pie_form = "*kʷoyséyeti"
    print(f"1. Proto-Indo-European Form:\t{pie_form}")
    print("   (Root *kʷoys- + Causative suffix -éye- + 3sg ending -ti)")
    print("-" * 70)

    # Step 2: Proto-Germanic (PGmc)
    # Grimm's Law: *kʷ > *hʷ, *t > *þ
    # Verner's Law: *s > *z (because PIE accent was on the suffix, not the root)
    # Vowel/Suffix change: *o > *a; *-éyeti becomes the weak verb suffix, resulting in *-iþ.
    # The full change is *kʷoyséyeti > *hʷaiziþ
    pgmc_form = "*hʷaiziþ"
    print(f"2. Proto-Germanic Changes:")
    print(f"   - Grimm's Law (kʷ > hʷ, t > þ) and Verner's Law (s > z) apply.")
    print(f"   - Vowel and suffix changes (o > a, -éyeti > -iþ).")
    print(f"   Proto-Germanic Form:\t\t{pgmc_form}")
    print("-" * 70)

    # Step 3: Old English (OE)
    # I-Umlaut: The *i in the suffix fronts the root vowel *ai > *ǣ.
    # Rhotacism: Intervocalic *z becomes *r.
    # Suffix leveling: The ending *-iþ becomes the standard Class I weak verb ending -eþ.
    oe_form = "*hwǣreþ"
    print(f"3. Old English Changes:")
    print(f"   - I-Umlaut (ai > ǣ, triggered by the following 'i').")
    print(f"   - Rhotacism (z > r).")
    print(f"   - Suffix regularization to '-eþ'.")
    print(f"   Old English Form:\t\t{oe_form}")
    print("-" * 70)
    
    # Step 4: Middle English (ME)
    # Vowel Shift: OE long *ǣ /æː/ becomes ME long open e /ɛː/, spelled <e>.
    # Spelling changes: *hw is now spelled <wh>, and *þ is spelled <th>.
    me_form = "whereth"
    print(f"4. Middle English Changes:")
    print(f"   - Vowel shift (ǣ > e).")
    print(f"   - Spelling changes (hw > wh, þ > th).")
    print(f"   Final Middle English Form:\t{me_form}")
    print("-" * 70)

derive_middle_english_form()