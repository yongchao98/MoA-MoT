import textwrap

def evaluate_mimicry_study_validity():
    """
    Analyzes the ecological validity of a proposed study on Bombus mimicry.
    The study uses human undergraduates to rate visual similarity.
    """

    # Helper function for formatted printing
    def print_section(title, content):
        print(f"\n--- {title.upper()} ---")
        # The textwrap.fill function wraps long lines of text
        print(textwrap.fill(content, width=80))

    # Step 1: Define the ecological function of the signal
    print_section(
        "Step 1: The Ecological Function of Mimicry",
        "Bumblebee (Bombus) mimicry is a form of Müllerian mimicry. Different, "
        "well-defended (stinging) species converge on a similar visual appearance. "
        "The function of this shared 'warning signal' is to more effectively "
        "teach predators to avoid them. The signal is therefore directed at predators."
    )

    # Step 2: Identify the signal receiver vs. the study's assessor
    print_section(
        "Step 2: The Signal Receiver vs. The Assessor",
        "The ecologically relevant 'signal receivers' are the natural predators of "
        "bumblebees, which are primarily birds. The proposed study, however, "
        "uses 'untrained undergraduates' (humans) as the assessors of similarity."
    )

    # Step 3: Compare the visual systems
    print_section(
        "Step 3: Comparing Visual Systems",
        "The validity of the study hinges on whether human vision is a good proxy "
        "for avian (bird) vision. There is a critical difference: \n\n"
        "1. Human Vision: Trichromatic. Humans have three types of color-sensitive "
        "cone cells (Red, Green, Blue) and cannot perceive ultraviolet (UV) light.\n\n"
        "2. Avian Vision: Tetrachromatic. Birds have four types of cone cells, "
        "one of which is sensitive to the UV part of the spectrum."
    )

    # Step 4: Analyze the implications of the visual system difference
    print_section(
        "Step 4: Implication of the Visual Difference",
        "Many insect color patterns, including the hairs on bumblebees, reflect UV "
        "light. These UV patterns are invisible to humans but are a salient and "
        "crucial part of the visual signal for birds. A study relying on human "
        "perception will completely miss this dimension of the mimicry signal. "
        "Furthermore, using static 2D images removes other potential cues predators "
        "might use, such as flight pattern and 3D shape."
    )

    # Step 5: Final Conclusion on Validity
    print("\n--- FINAL CONCLUSION ---")
    print(
        "\nIs the approach valid? \n\n"
        "No. The approach is not ecologically valid.\n"
        "Because human vision is a poor proxy for the vision of the predators that "
        "drive the evolution of mimicry, the resulting similarity clusters would be "
        "based on incomplete and potentially misleading information. The study would "
        "measure human-perceived similarity, not the predator-perceived similarity "
        "that actually matters for the mimicry syndrome's function."
    )


if __name__ == "__main__":
    evaluate_mimicry_study_validity()