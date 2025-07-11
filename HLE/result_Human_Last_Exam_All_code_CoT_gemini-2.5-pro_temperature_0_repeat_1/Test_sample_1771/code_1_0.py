def trace_etymology():
    """
    This function prints the etymological derivation of a PIE word
    through Germanic languages.
    """
    pie_form = "*seh₂gieti"
    pgmc_form = "*sōkiþ"
    on_form = "søkir"
    oe_form = "sēċeþ"
    ohg_form = "suohhit"

    print(f"The derivation from PIE {pie_form} '(s)he seeks/gives a sign' is as follows:\n")

    # Print each stage of the derivation
    print(f"1. Proto-Indo-European: {pie_form}")
    print(f"2. Proto-Germanic:      {pgmc_form}")
    print(f"3. Old Norse:           {on_form}")
    print(f"4. Old English:         {oe_form}")
    print(f"5. Old High German:     {ohg_form}\n")

    # Print the final derivation chain as a single equation
    print("The complete derivation chain is:")
    final_equation = f"{pie_form} > {pgmc_form} > (ON {on_form}, OE {oe_form}, OHG {ohg_form})"
    
    # As requested, output each part of the final equation
    print(f"PIE form: {pie_form}")
    print(f"Proto-Germanic form: {pgmc_form}")
    print(f"Old Norse form: {on_form}")
    print(f"Old English form: {oe_form}")
    print(f"Old High German form: {ohg_form}")


trace_etymology()