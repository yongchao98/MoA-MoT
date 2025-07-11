def solve_lebesgue_integrability():
    """
    Analyzes a list of function types and determines which are necessarily
    Lebesgue integrable.

    The definition of a Lebesgue integrable function 'f' is:
    1. f is measurable.
    2. The integral of its absolute value is finite (∫|f|dμ < ∞).

    This function encodes the reasoning for each choice and prints the final
    answer as a single string.
    """

    choices = {
        'A': "A bounded function: No. May not be measurable.",
        'B': "A bounded measurable function: No. Domain may have infinite measure.",
        'C': "A measurable function: No. Integral of absolute value may be infinite.",
        'D': "A continuous function: No. Domain may have infinite measure.",
        'E': "A measurable function on [a,b]: No. Integral of absolute value may be infinite (e.g., 1/x on (0,1]).",
        'F': "A continuous function on [a,b]: Yes. A continuous function on a compact set is bounded, so ∫|f| ≤ M(b-a) < ∞.",
        'G': "A bounded function on [a,b]: No. May not be measurable.",
        'H1': "A function whose absolute value is integrable: No. The function f itself may not be measurable.",
        'I': "A function whose absolute value is integrable on [a,b]: No. The function f itself may not be measurable.",
        'J': "A continuous function on (a,b): No. Integral may diverge at an endpoint (e.g., 1/(x-a)).",
        'H2': "A bounded function on (a,b): No. May not be measurable.",
        'K': "A measurable function on (a,b): No. Integral may diverge (e.g., 1/(x-a)).",
        'L': "A measurable function whose absolute value is integrable: Yes. This is the definition of Lebesgue integrability.",
        'M': "A bounded continuous function on (a,b): Yes. f is measurable and ∫|f| ≤ M(b-a) < ∞."
    }

    # The actual letters from the prompt. Note 'H' appears twice with different meanings.
    # We will map our analysis keys to the prompt's letters.
    prompt_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'H', 'K', 'L', 'M']
    
    correct_choices = []
    
    # Logic to find the correct choices based on the analysis
    if "Yes" in choices['F']:
        correct_choices.append('F')
    if "Yes" in choices['L']:
        correct_choices.append('L')
    if "Yes" in choices['M']:
        correct_choices.append('M')

    # The problem asks for the answer as a string of letters in order.
    # F, L, M are already in alphabetical/list order.
    answer = "".join(sorted(correct_choices))
    
    print(answer)

solve_lebesgue_integrability()