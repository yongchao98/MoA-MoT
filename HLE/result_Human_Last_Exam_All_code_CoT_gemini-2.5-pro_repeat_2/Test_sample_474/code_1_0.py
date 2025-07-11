def evaluate_research_design():
    """
    This function evaluates the validity of a research design for studying
    bumblebee mimicry syndromes using human perception.
    """

    # --- Step 1: Define the core components of the experiment ---
    undergraduate_count = 20
    print(f"The experiment uses {undergraduate_count} untrained human undergraduates to rank bee similarity.")
    print("The goal is to cluster bees into mimicry syndromes based on these rankings.")
    print("-" * 50)

    # --- Step 2: Define the ecological function of mimicry ---
    print("ASSESSING THE ECOLOGICAL FUNCTION:")
    print("The ecological function of Müllerian mimicry in bumblebees is to present a memorable, common warning signal to potential predators.")
    print("The effectiveness of the signal depends entirely on the perception of the PREDATOR, not humans.")
    print("-" * 50)

    # --- Step 3: Compare the visual systems of the observers ---
    print("COMPARING OBSERVERS: Human vs. Natural Predator")
    print("\nHuman Observers (the undergraduates):")
    print(" - Vision: Trichromatic (sensitive to Red, Green, Blue light).")
    print(" - Cannot perceive ultraviolet (UV) light.")

    print("\nNatural Predators (e.g., Birds):")
    print(" - Vision: Tetrachromatic (sensitive to Red, Green, Blue, AND Ultraviolet light).")
    print(" - Many floral patterns and insect colorations have UV components invisible to humans.")
    print("-" * 50)

    # --- Step 4: Draw the conclusion ---
    print("CONCLUSION ON VALIDITY:")
    print("The central flaw in the research design is using an observer (human) whose sensory perception is fundamentally different from the ecologically relevant observer (predator).")
    print("\nTwo bee species that appear identical to the {undergraduate_count} humans might look completely different to a bird due to differences in their UV coloration.".format(undergraduate_count=undergraduate_count))
    print("Conversely, two species that humans rank as different might be functionally identical as a warning signal to a bird.")
    print("\nTherefore, the data on human-perceived similarity is not a valid proxy for the functional similarity of the mimicry signal in its natural, ecological context.")
    print("-" * 50)

if __name__ == '__main__':
    evaluate_research_design()
    print("\nFINAL ANSWER: Is the approach valid considering the ecological function?")
    print("No, the approach is not valid.")
