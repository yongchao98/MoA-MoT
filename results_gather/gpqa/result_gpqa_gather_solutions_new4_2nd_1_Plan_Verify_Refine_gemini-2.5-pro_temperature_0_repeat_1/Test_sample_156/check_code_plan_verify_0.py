import re

def check_answer_correctness(question, llm_answer_text):
    """
    Checks the correctness of the LLM's answer by verifying it against the question's constraints.

    Args:
        question (str): The original question text.
        llm_answer_text (str): The full text of the LLM's answer, including the final choice.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """

    # --- Step 1: Define the constraints from the question ---
    constraints = {
        "kit_type": "molecular",
        "pathogen_type": "retrovirus",
        "detection_speed": "quick",
        "purpose": "outbreak"
    }

    # --- Step 2: Define the properties of each option based on scientific principles ---
    # These properties will be used to evaluate if an option meets the constraints.
    options_analysis = {
        'A': {
            "description": "Identify by symptoms, then design a nested PCR kit.",
            "is_molecular": True,
            "handles_retrovirus": True, # PCR is a valid technique for retroviral cDNA
            "is_quick_for_outbreak": True, # PCR is fast
            "is_scientifically_valid_method": False, # Cannot design a PCR test from symptoms
            "reason_for_failure": "The identification method (using symptoms) is scientifically invalid for designing a sequence-specific molecular test."
        },
        'B': {
            "description": "Identify IgG antibodies, then develop an ELISA kit.",
            "is_molecular": False, # ELISA is serological/immunological
            "handles_retrovirus": False, # Does not detect genetic material
            "is_quick_for_outbreak": False, # IgG is a late-stage marker, not for early/quick detection
            "is_scientifically_valid_method": True, # As a serological test, it's valid, but it doesn't fit the question
            "reason_for_failure": "The proposed kit is not a 'molecular' kit, and targeting IgG is not suitable for 'quick' early detection during an outbreak."
        },
        'C': {
            "description": "Identify by cDNA sequencing, then develop a real-time PCR kit.",
            "is_molecular": True,
            "handles_retrovirus": True, # cDNA sequencing correctly handles the RNA genome
            "is_quick_for_outbreak": True, # Real-time PCR is the gold standard for speed and accuracy
            "is_scientifically_valid_method": True,
            "reason_for_failure": None # This is the correct option
        },
        'D': {
            "description": "Identify by DNA sequencing, then develop a PCR kit.",
            "is_molecular": True,
            "handles_retrovirus": False, # Incorrectly assumes a DNA genome for a retrovirus
            "is_quick_for_outbreak": True,
            "is_scientifically_valid_method": False, # The initial step is biologically incorrect for a retrovirus
            "reason_for_failure": "The method fails to account for the fact that a retrovirus has an RNA genome, not a DNA genome. It misses the essential reverse transcription (RNA to cDNA) step."
        }
    }

    # --- Step 3: Identify the correct option based on the analysis ---
    correct_option = None
    for option, properties in options_analysis.items():
        if properties["reason_for_failure"] is None:
            correct_option = option
            break

    # --- Step 4: Extract the LLM's chosen answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: Could not find the final answer in the format <<<A>>>, <<<B>>>, etc."
    
    chosen_option = match.group(1)

    # --- Step 5: Compare the chosen answer with the correct answer and provide feedback ---
    if chosen_option == correct_option:
        return "Correct"
    else:
        failure_reason = options_analysis[chosen_option]["reason_for_failure"]
        return (f"Incorrect. The final answer was {chosen_option}, but the correct answer is {correct_option}. "
                f"The chosen option is wrong because: {failure_reason}")


# Provided question and LLM answer
question_text = """
There has been an outbreak of an viral infectious disease in your city. You have been asked to design a molecular diagnostic kit for quick detection of this retrovirus. How would you go about designing this kit?

A) First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis.
B) First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick 
 and accurate diagnosis.
C) First you identify the virus by performing cDNA sequencing.  You then develop a real time PCR kit for quick and accurate diagnosis.
D) First you identify the virus by performing DNA sequencing. You then develop a PCR kit for  quick and accurate diagnosis
"""

llm_final_answer = """
Here is a step-by-step analysis of the question, incorporating all the careful points derived from the candidate answers.

### **1. Deconstruct the Core Task**

The question asks for the design of a **molecular diagnostic kit** for a **retrovirus** outbreak, emphasizing the need for **quick and accurate detection**.

*   **Careful Point 1: "Molecular Diagnostic Kit"**: This is a critical constraint. A molecular test must detect the pathogen's genetic material (nucleic acids, i.e., RNA or DNA). This is fundamentally different from serological or immunological tests, which detect the host's immune response (e.g., antibodies) or viral proteins (antigens).

*   **Careful Point 2: "Retrovirus"**: This is the key biological detail. A retrovirus has an RNA genome. A defining feature is its use of the enzyme reverse transcriptase to create a DNA copy of its RNA genome (a process called reverse transcription). The resulting DNA is called complementary DNA (cDNA). Any molecular method based on DNA amplification or sequencing must first account for this RNA-to-cDNA conversion step.

*   **Careful Point 3: "Quick and Accurate Detection"**: This requirement is crucial for managing an outbreak. "Quick" implies a rapid turnaround time for results. "Accurate" means the test must be highly sensitive (correctly identifying those with the virus, even with low viral loads) and specific (correctly identifying those without the virus). This points toward a method suitable for early detection.

### **2. Evaluate Each Option Against the Core Task**

**A) First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis.**
*   **Careful Point:** The identification method is scientifically invalid. Clinical symptoms like fever or cough are non-specific and can be caused by numerous different pathogens. It is impossible to determine the precise genetic sequence of a virus—which is required to design specific PCR primers—based on symptoms alone.

**B) First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick and accurate diagnosis.**
*   **Careful Point 1:** This describes a serological (immunological) test, not a **molecular** one. It fails to meet the primary requirement of the task.
*   **Careful Point 2:** This method is not suitable for **quick**, early detection. IgG antibodies are part of the long-term immune response and typically take several days to weeks to become detectable after infection. This creates a "window period" where an infected person is contagious but would test negative, making it a poor choice for managing an acute outbreak.

**C) First you identify the virus by performing cDNA sequencing. You then develop a real time PCR kit for quick and accurate diagnosis.**
*   **Careful Point 1:** The identification step, **"cDNA sequencing"**, is correct. It acknowledges that the retrovirus has an RNA genome that must first be converted to complementary DNA (cDNA) before it can be sequenced. This is the definitive way to identify the virus's unique genetic code.
*   **Careful Point 2:** The detection method, a **"real time PCR kit"**, is the modern gold standard. Real-time PCR (specifically, RT-qPCR, which incorporates the reverse transcription step) is extremely fast, highly sensitive, and specific. It perfectly fits the "quick and accurate" requirement for an outbreak scenario. This option describes a complete, scientifically sound, and state-of-the-art workflow.

**D) First you identify the virus by performing DNA sequencing. You then develop a PCR kit for quick and accurate diagnosis.**
*   **Careful Point:** This option is based on a fundamental biological error. It assumes the retrovirus has a DNA genome. Attempting to perform direct **"DNA sequencing"** on the virus particle would fail, as its genome is RNA. This option misses the essential reverse transcription step.

### **Conclusion**

Based on the comprehensive analysis, Option C is the only one that presents a scientifically correct, logical, and modern workflow that satisfies all conditions of the question. It correctly identifies the nature of a retrovirus, the necessary steps to analyze its genetic material, and the most appropriate technology for rapid and accurate diagnosis.

<<<C>>>
"""

# Run the check
result = check_answer_correctness(question_text, llm_final_answer)
print(result)