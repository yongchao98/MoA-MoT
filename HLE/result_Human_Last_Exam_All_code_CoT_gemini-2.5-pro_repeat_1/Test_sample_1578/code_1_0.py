import textwrap

def solve_mirror_puzzle():
    """
    This script analyzes the Japanese puzzle of the Mirror and the Oni
    to determine the best method for detecting the demon.
    """

    print("--- The Mirror and the Oni Puzzle ---")
    print("The task is to distinguish a true reflection from a demon's illusion.\n")

    # The chosen method is 'F'. We will analyze it in detail.
    chosen_method = 'F'
    method_text = "Use a second mirror to reflect the alleged first mirror to catch the demon."

    print(f"Chosen Method for Analysis: {chosen_method}")
    print(textwrap.fill(f"Description: {method_text}", width=70))
    print("-" * 35)

    print("\nLogical Breakdown of the Solution:\n")

    # Step 1: The setup
    print("1. The Scenario:")
    print("   - You are looking into the primary mirror (Mirror 1).")
    print("   - A figure that looks exactly like you copies your every move.")
    print("   - This figure is either your true reflection or a demon mimicking you.\n")

    # Step 2: The Test
    print("2. The Test with a Second Mirror (Mirror 2):")
    print("   - You hold a second, smaller mirror at an angle.")
    print("   - You position Mirror 2 so that it shows you the reflection from Mirror 1.\n")

    # Step 3: Expected Outcomes
    print("3. Two Possible Outcomes:\n")

    # Case A: It's a true reflection
    print("   Case A: If it is a REAL REFLECTION...")
    explanation_a = (
        "The laws of optics apply. Light from the back of your head hits Mirror 1, "
        "bounces to Mirror 2, and then reflects to your eyes. In Mirror 2, "
        "you will see the reflection of the reflection of the back of your own head."
    )
    print(textwrap.fill(explanation_a, width=65, initial_indent='     ', subsequent_indent='     '))

    # Case B: It's a demon
    print("\n   Case B: If it is a DEMON...")
    explanation_b = (
        "The demon is a physical being standing where your reflection should be, facing you. "
        "It is not a product of optics. When you look into Mirror 2, it will reflect "
        "the image of what is in Mirror 1, which is the demon. Since the demon is facing you, "
        "you will see the demon's face in Mirror 2."
    )
    print(textwrap.fill(explanation_b, width=65, initial_indent='     ', subsequent_indent='     '))

    # Step 4: The Conclusion
    print("\n4. Conclusion:")
    conclusion = (
        "The two outcomes are different. The demon cannot replicate the correct "
        "perspective-based reflection from a different angle while simultaneously "
        "facing you to maintain the primary illusion. This geometric paradox "
        "exposes the demon."
    )
    print(textwrap.fill(conclusion, width=65, initial_indent='   '))
    print("\nTherefore, method 'F' is the most effective solution.")


# Execute the solver
solve_mirror_puzzle()

# There is no equation in this logic puzzle. The instruction to output numbers
# from an equation does not apply here. The final answer is the letter 'F'.
