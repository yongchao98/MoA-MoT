import textwrap

def explain_poker_gto_concepts():
    """
    Analyzes subtle betting reasons in poker under GTO assumptions.
    """

    print("Analyzing subtle betting reasons under Game-Theoretically Optimal (GTO) assumptions:")
    print("The core assumption is that both you and your opponent play a perfect, unexploitable strategy.")
    print("We will evaluate reasons 3, 4, and 5.\n")

    reasons = {
        3: {
            "name": "Denying equity to drawing hands.",
            "disappears": False,
            "explanation": (
                "This reason remains valid. Denying equity is a fundamental part of GTO. "
                "A 'protection bet' forces hands with significant drawing equity to either fold or pay an "
                "unfavorable price to continue. This is a mathematically sound component of a "
                "balanced betting strategy and does NOT disappear."
            )
        },
        4: {
            "name": "Gaining information about your opponent's hand.",
            "disappears": True,
            "explanation": (
                "This reason disappears. A GTO opponent's response to a bet is perfectly balanced. "
                "Their calling and raising ranges contain a calculated mix of value hands and bluffs. "
                "Therefore, their action does not give away any 'exploitable' information. "
                "Betting to 'see where you're at' is an exploitative concept that is "
                "ineffective against a perfect GTO player."
            )
        },
        5: {
            "name": "Avoiding revealing your own hand in a showdown.",
            "disappears": True,
            "explanation": (
                "This reason disappears. A GTO opponent cannot exploit the knowledge of your single, "
                "specific hand revealed at showdown. Your strategy is defined by balanced ranges "
                "for every action. A GTO opponent already knows these ranges and cannot profitably "
                "alter their own perfect strategy based on seeing one hand from that range. "
                "The fear of giving away information about your tendencies is irrelevant."
            )
        }
    }

    disappearing_reasons = []
    for num, details in reasons.items():
        print(f"--- Analysis of Reason {num}: {details['name']} ---")
        if details["disappears"]:
            status = "DISAPPEARS"
            disappearing_reasons.append(str(num))
        else:
            status = "DOES NOT disappear"
        
        print(f"Status: {status}")
        print("Explanation:")
        # textwrap makes the output cleaner
        for line in textwrap.wrap(details["explanation"], width=70):
            print(f"  {line}")
        print()

    print("="*40)
    final_conclusion = " and ".join(disappearing_reasons)
    print(f"Conclusion: The reasons that disappear under GTO are {final_conclusion}.")
    print("="*40)


if __name__ == "__main__":
    explain_poker_gto_concepts()