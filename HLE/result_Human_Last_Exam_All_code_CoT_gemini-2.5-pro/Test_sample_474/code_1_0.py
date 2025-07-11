import textwrap

def evaluate_research_method():
    """
    Analyzes the validity of a proposed research method for clustering
    Bombus species by mimicry syndrome.
    """

    print("Analysis of Research Method: Clustering Bombus by Mimicry Syndrome")
    print("=" * 70)
    
    # Define the question
    question = "Is it a valid approach to use untrained undergraduates to visually rank Bombus species similarity from field images to cluster them by mimicry syndrome?"
    print(f"Question: {textwrap.fill(question, 70)}\n")

    # Step 1: Explain the ecological function
    print("Step 1: The Ecological Function of Bumblebee Mimicry")
    print("-" * 70)
    explanation1 = "Bumblebee mimicry is a form of Müllerian mimicry. Multiple well-defended species (they can sting) evolve to share a similar warning color pattern."
    explanation2 = "The function of this shared signal is to reduce the per-species cost of teaching predators to avoid them. The key signal 'receiver' in this ecological context is the PREDATOR (e.g., birds)."
    print(textwrap.fill(explanation1, 70))
    print(textwrap.fill(explanation2, 70))
    print("\n")

    # Step 2: Evaluate the proposed methodology against the ecological function
    print("Step 2: Evaluating the Proposed Method's Validity")
    print("-" * 70)
    
    flaw1_title = "Flaw 1: The Observer is Ecologically Irrelevant."
    flaw1_desc = "The study uses humans, but the mimicry signal is directed at natural predators. Any resulting clusters would be based on human perception, not predator perception, making them irrelevant to the ecological function of mimicry."
    print(f"-> {flaw1_title}")
    print(f"   {textwrap.fill(flaw1_desc, 67, initial_indent=' ')}")
    print()

    flaw2_title = "Flaw 2: Human and Predator Vision Are Fundamentally Different."
    flaw2_desc = ("Humans have trichromatic vision (Red, Green, Blue). Key bumblebee predators, like birds, have tetrachromatic vision, allowing them to see in the ultraviolet (UV) spectrum. Bumblebee color patterns often have significant UV components that are completely invisible to humans. Two species that look similar to a human may look very different to a bird, invalidating the entire approach.")
    print(f"-> {flaw2_title}")
    print(f"   {textwrap.fill(flaw2_desc, 67, initial_indent=' ')}")
    print()

    flaw3_title = "Flaw 3: 'Untrained' Observers Introduce Uncontrolled Bias."
    flaw3_desc = "Untrained individuals may base their similarity rankings on features unrelated to the aposematic signal, such as body size, posture, image quality, or even the flower the bee is on, introducing significant noise and bias."
    print(f"-> {flaw3_title}")
    print(f"   {textwrap.fill(flaw3_desc, 67, initial_indent=' ')}")
    print("\n")

    # Step 3: Final Conclusion
    print("Step 3: Conclusion")
    print("-" * 70)
    conclusion = "The proposed approach is NOT valid. It fails to account for the most critical component of the mimicry system: the visual perception of the natural predators that the signals are intended for. The results would not reflect the true biological mimicry rings."
    print(textwrap.fill(conclusion, 70))

if __name__ == "__main__":
    evaluate_research_method()