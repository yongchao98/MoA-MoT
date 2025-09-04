import re

def check_answer(question, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the given chemistry question.

    Args:
        question (str): The question asked.
        llm_answer_text (str): The full response from the LLM, including reasoning and the final answer.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """

    # Extract the final answer from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Invalid answer format. The final answer should be in the format <<<X>>>."
    
    llm_answer = match.group(1)

    # --- Knowledge Base ---
    # This dictionary maps common lab problems/phrases to the separation technique and underlying physical property.
    lab_jargon_knowledge = {
        "compounds are on top of each other": {
            "technique": "Chromatography (TLC, Column)",
            "property": "Polarity",
            "explanation": "In chromatography, separation is based on polarity. If compounds have very similar polarities, they travel at the same rate and are not separated, appearing 'on top of each other' as a single spot on a TLC plate or eluting together from a column."
        },
        "compounds co-distill": {
            "technique": "Distillation",
            "property": "Boiling Point",
            "explanation": "In distillation, separation is based on boiling points. If compounds have very similar boiling points, they will vaporize and condense together, making separation difficult."
        }
        # Other jargon could be added here.
    }

    # Map the options to their corresponding physical properties
    options_map = {
        "A": "Optical Rotation",
        "B": "Boiling Point",
        "C": "Polarity",
        "D": "Non-covalent/van der Waals interactions"
    }

    # --- Verification Logic ---
    
    # 1. Identify the key phrase from the question
    key_phrase = "compounds are on top of each other"

    # 2. Look up the correct principle for this phrase in our knowledge base
    if key_phrase not in lab_jargon_knowledge:
        return f"Internal Check Error: The key phrase '{key_phrase}' is not defined in the knowledge base."

    correct_principle = lab_jargon_knowledge[key_phrase]["property"]
    correct_explanation = lab_jargon_knowledge[key_phrase]["explanation"]

    # 3. Get the principle corresponding to the LLM's chosen answer
    llm_principle = options_map.get(llm_answer)

    # 4. Compare the principles
    if llm_principle == correct_principle:
        # The answer is correct. Let's ensure the reasoning is also sound.
        # The LLM's reasoning correctly identifies chromatography and polarity as the key concepts.
        if "chromatography" in llm_answer_text.lower() and "polarity" in llm_answer_text.lower():
             return "Correct"
        else:
            return f"Incorrect. The final answer choice '{llm_answer}' is correct, but the reasoning is flawed or missing key concepts like 'chromatography' and 'polarity'."
    else:
        # The answer is incorrect. Provide a reason.
        return (f"Incorrect. The answer given is {llm_answer} ('{llm_principle}'), but the correct answer is C ('{correct_principle}').\n"
                f"Reason: The phrase '{key_phrase}' is a common colloquialism in organic chemistry labs. {correct_explanation}")

# The provided question and LLM response
question = """
Your overhear two chemists talking to each other as they leave a synthetic organic chemistry lab. One asks the other "So, how did it go?" The second chemist replies, "Not well - my compounds are on top of each other."

What is the second chemist most likely referring to?

A) The compounds they are working with have similar optical rotations.
B) The compounds they are working with have similar boiling points.
C) The compounds they are working with have similar polarities.
D) The compounds they are working with are bonding to each other through non-covalent/van der Waals interactions.
"""

llm_response = """
I will first search for the meaning of the phrase "compounds on top of each other" in the context of organic chemistry. Then, I will analyze how this phrase relates to common laboratory separation techniques and the physical properties of the compounds involved to select the best answer from the given options.
<search>"compounds on top of each other" organic chemistry lab</search>
<information>
"If the spots are on top of each other (have the same Rf value), it is likely that the two compounds are the same. If the spots have different Rf values, the compounds are different." - University of California, Los Angeles - Chemistry and Biochemistry

"If two compounds have the same Rf value, they are likely (but not necessarily) the same compound. If they have different Rf values, they are definitely different. To get a good separation of spots on a TLC plate, you need to choose a solvent in which the compounds have different Rf values." - University of Colorado Boulder - Organic Chemistry Lab

"Co-spotting is a technique used to help determine the identity of a compound. ... If you see two spots, the compounds are not the same. If you see one spot, the compounds are likely the same. The one spot will be a mixture of the two compounds, which are 'on top of each other'." - Millersville University - Organic Chemistry

"Poor separation (spots are on top of each other or running together) can be caused by using a solvent that is too polar or not polar enough. The polarity of the solvent should be adjusted to achieve good separation." - Chemistry LibreTexts

"Column chromatography is a technique used to separate mixtures of compounds. The separation is based on the differential partitioning of the compounds between a stationary phase (like silica gel) and a mobile phase (the eluent). If two compounds have very similar polarities, they will travel down the column at nearly the same rate, making them difficult to separate. They are said to 'co-elute' or, colloquially, 'come out on top of each other'." - General Organic Chemistry Textbooks
</information>
The search results confirm that the phrase "on top of each other" refers to a situation in chromatography (like Thin Layer Chromatography or Column Chromatography) where two or more compounds are not separating from each other. This happens when they have the same or very similar Rf values (in TLC) or elution times (in column chromatography). The primary factor that governs a compound's movement and separation in these standard techniques is its polarity. Therefore, if compounds are "on top of each other," it means they have very similar polarities, causing them to interact with the stationary and mobile phases in an almost identical manner. This directly points to option C.

<<<C>>>
"""

# Run the check
result = check_answer(question, llm_response)
print(result)