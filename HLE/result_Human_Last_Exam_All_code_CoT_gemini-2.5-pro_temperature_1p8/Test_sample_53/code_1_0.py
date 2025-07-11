def derive_middle_english_form():
    """
    This function traces a PIE root through its historical sound changes
    to its theoretical form in Middle English.
    """
    # 1. Start with the PIE root and apply o-grade.
    # The root is *kʷeys-, o-grade is *kʷoys-.
    # The causative form *kʷoys-éye- is the basis for Verner's Law.
    pie_stem = "kʷoys"
    print(f"1. PIE o-grade causative stem: *{pie_stem}-")

    # 2. Apply sound changes for Proto-Germanic.
    # Grimm's Law: kʷ > hw, o > a.
    # Verner's Law: s > z (due to accent on suffix).
    pgmc_stem = "hw" + pie_stem[2:]  # kʷ > hw
    pgmc_stem = pgmc_stem.replace('o', 'a')
    pgmc_stem = pgmc_stem.replace('s', 'z') # Verner's Law
    print(f"2. Proto-Germanic stem (after Grimm's & Verner's Laws): *{pgmc_stem}-")

    # 3. Apply sound changes for Old English.
    # West Germanic Rhotacism: z > r.
    # I-Umlaut: ai > ǣ.
    # Add 3rd person singular present suffix: -eþ.
    oe_stem = pgmc_stem.replace('z', 'r') # Rhotacism
    oe_stem = oe_stem.replace('ai', 'ǣ')  # I-Umlaut
    oe_verb_form = oe_stem + "eþ"
    print(f"3. Old English 3rd person singular present form: {oe_verb_form}")

    # 4. Apply sound changes for Middle English.
    # Spelling: hw- > wh-.
    # Vowel shift: ǣ > e.
    # Suffix remains as -eth.
    me_verb_form = "wh" + oe_verb_form[2:] # Spelling change
    me_verb_form = me_verb_form.replace('ǣ', 'e')
    # The suffix -eþ becomes -eth. This replacement is safe here.
    me_verb_form = me_verb_form.replace('eþ', 'eth')
    print(f"4. Middle English form: {me_verb_form}")

    print(f"\nTherefore, the third person singular present verbal form, meaning 'he shows', would be '{me_verb_form}'.")
    return me_verb_form

# Execute the derivation and capture the final answer.
final_answer = derive_middle_english_form()

# The final answer in the requested format.
print(f"\n<<<he {final_answer}>>>")