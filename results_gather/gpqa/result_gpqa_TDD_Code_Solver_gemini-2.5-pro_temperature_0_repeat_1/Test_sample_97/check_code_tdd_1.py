import re

def check_chemistry_answer(llm_answer: str):
    """
    Checks the correctness of a given answer to a multiple-choice chemistry question.

    The question is:
    Identify the starting material, A, in the following reaction.
    A + a methyleneruthenium compound + 1-propene ---> 1-(prop-1-en-1-yl)-2-vinylcyclopentane

    A) 2-methylbicyclo[3.1.0]hex-2-ene
    B) 2-methyl-3-methylenebicyclo[2.1.0]pentane
    C) bicyclo[3.2.0]hept-6-ene
    D) 1,2-dimethylenecyclopentane

    This function determines the correct answer based on established chemical reaction
    mechanisms and then compares it to the provided answer.
    """

    # Step 1: Determine the correct answer based on chemical principles.
    # The reaction described is a tandem Ring-Opening Metathesis / Cross-Metathesis (ROCM).
    # This reaction type requires a strained cyclic alkene as a starting material that, upon
    # opening, can lead to the product's carbon skeleton.
    #
    # - The product is 1-(prop-1-en-1-yl)-2-vinylcyclopentane. It has a cyclopentane core
    #   with two adjacent substituents: a vinyl group (-CH=CH2) and a propenyl group (-CH=CH-CH3).
    # - The reaction involves a ruthenium catalyst and 1-propene, which are characteristic
    #   reagents for olefin metathesis.
    #
    # - Analysis of Option C, bicyclo[3.2.0]hept-6-ene:
    #   1. Ring-Opening Metathesis (ROM): The strained cyclobutene ring of bicyclo[3.2.0]hept-6-ene
    #      is opened by the ruthenium catalyst. This creates a cis-1,2-disubstituted cyclopentane
    #      intermediate which is itself a ruthenium carbene: [Ru]=CH-(cyclopentane)-CH=CH2.
    #   2. Cross-Metathesis (CM): This intermediate carbene reacts with 1-propene (CH3-CH=CH2).
    #      The [Ru]=CH- group of the carbene is exchanged for the CH3-CH= group from propene.
    #   3. The final product is CH3-CH=CH-(cyclopentane)-CH=CH2, which is
    #      1-(prop-1-en-1-yl)-2-vinylcyclopentane. This pathway perfectly matches the given transformation.
    #
    # - Other options are incorrect as their ring-opening or cross-metathesis reactions
    #   would not yield the specific product structure.
    #
    # Therefore, the correct option is 'C'.
    correct_option = 'C'
    options = ['A', 'B', 'C', 'D']

    # Step 2: Parse the provided LLM answer to find the selected option.
    # The provided answer is conversational text and does not select an option.
    # We will check for a clear answer format like <<<C>>> or a standalone letter.
    
    selected_option = None
    # A robust regex to find a single capital letter from the options, not part of a word.
    match = re.search(r'(?< ![a-zA-Z]) ([A-D])(? ![a-zA-Z]) ', llm_answer)
    if match:
        selected_option = match.group(1)
    
    # Also check for the specific format <<<X>>>
    match_specific = re.search(r'<<<([A-D])>>>', llm_answer)
    if match_specific:
        selected_option = match_specific.group(1)

    # Step 3: Compare the found option with the correct option and return the result.
    if selected_option is None:
        return "Incorrect. The provided response is not a valid answer. It is conversational text that fails to select one of the options (A, B, C, or D)."

    if selected_option == correct_option:
        return "Correct"
    else:
        return f"Incorrect. The selected option is {selected_option}, but the correct option is {correct_option}. The reaction is a tandem Ring-Opening-Cross-Metathesis (ROCM) starting from bicyclo[3.2.0]hept-6-ene."

# # Example of how to use the checker with the provided LLM answer.
# llm_response_text = """Excellent! The test-driven development process successfully identified the correct starting material, and the checker confirmed the result.
#
# I am ready for the next question.
# Your checking code:"""
#
# result = check_chemistry_answer(llm_response_text)
# print(result)
# # Expected output: "Incorrect. The provided response is not a valid answer. It is conversational text that fails to select one of the options (A, B, C, or D)."
