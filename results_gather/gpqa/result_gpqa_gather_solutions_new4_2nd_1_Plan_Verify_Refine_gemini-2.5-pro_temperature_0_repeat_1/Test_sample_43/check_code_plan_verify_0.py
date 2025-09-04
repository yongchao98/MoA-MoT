import re

def check_answer(question_text, llm_answer_text):
    """
    Checks the correctness of the LLM's answer to a multiple-choice question about mitochondrial assays.

    The function encodes the biological principles behind each experimental option to determine which one is not helpful for studying mitochondrial function.

    Args:
        question_text (str): The original question text.
        llm_answer_text (str): The full text of the LLM's response, including the final answer in <<<>>> format.

    Returns:
        str: "Correct" if the answer is correct, or a string explaining the error if it's incorrect.
    """

    # --- Step 1: Define the biological validity of each option ---
    # The question asks for the experiment that will NOT help.
    # We will define the validity of each experiment. The correct answer is the one that is NOT valid.
    
    experimental_validity = {
        'A': {
            'valid': True,
            'name': "Confocal fluorescence microscopy after Mito-RTP staining",
            'reason': "This is a valid method. 'Mito-' probes are designed to measure specific parameters like temperature or ROS directly within the mitochondria, which are direct indicators of metabolic activity and stress."
        },
        'B': {
            'valid': True,
            'name': "Transformation of cells with recombinant luciferase and luminometer reading",
            'reason': "This is a valid method. It measures ATP levels, which are the primary output of mitochondrial energy production. Changes in ATP reflect changes in mitochondrial function."
        },
        'C': {
            'valid': True,
            'name': "Flow cytometry after labeling with JC-1 dye",
            'reason': "This is a valid method. JC-1 dye measures the mitochondrial membrane potential, a critical indicator of mitochondrial health and the integrity of the electron transport chain."
        },
        'D': {
            'valid': False,
            'name': "Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit",
            'reason': "This method is biologically flawed. Glucose uptake occurs at the cell's plasma membrane, not the mitochondrial membrane. Mitochondria import pyruvate (the product of glycolysis in the cytoplasm), not glucose itself. Therefore, this assay is nonsensical for isolated mitochondria."
        }
    }

    # --- Step 2: Identify the correct answer based on the defined logic ---
    correct_option = None
    for option, details in experimental_validity.items():
        if not details['valid']:
            correct_option = option
            break
    
    if correct_option is None:
        return "Error in checker logic: No single invalid option was identified."

    # --- Step 3: Extract the LLM's chosen answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Could not find the final answer in the required format '<<<A>>>', '<<<B>>>', etc."
    
    llm_option = match.group(1)

    # --- Step 4: Compare the LLM's answer with the correct answer and provide feedback ---
    if llm_option == correct_option:
        return "Correct"
    else:
        # The LLM chose a valid experiment as the "unhelpful" one.
        llm_choice_details = experimental_validity[llm_option]
        correct_choice_details = experimental_validity[correct_option]
        
        reasoning = (
            f"Incorrect. The chosen answer '{llm_option}' is wrong.\n"
            f"The question asks for the experiment that will NOT help. "
            f"Option {llm_option} ({llm_choice_details['name']}) is actually a helpful experiment. {llm_choice_details['reason']}\n"
            f"The correct answer is '{correct_option}' ({correct_choice_details['name']}). {correct_choice_details['reason']}"
        )
        return reasoning

# Provided question and answer text from the user prompt
question = """
A eukaryotic cell evolved a mechanism to turn macromolecular building blocks into energy. The process occurs in mitochondria, which are cellular energy factories. In the series of redox reactions, the energy from food is stored between the phosphate groups and used as a universal cellular currency. The energy-laden molecules are shuttled out of the mitochondrion to serve in all cellular processes. You discovered a new anti-diabetes drug and want to investigate whether it has an effect on the mitochondria. You set up a range of experiments with your HEK293 cell line. Which of the experiments listed below will not help you discover the mitochondrial role of your drug:

A) Confocal fluorescence microscopy after Mito-RTP staining of the cells
B) Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant
C) Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide
D) Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit
"""

llm_answer = """
Here is a step-by-step analysis of the question to determine the final answer.

### Step 1: Understand the Goal
The objective is to identify which of the four listed experiments is **not** suitable for investigating a drug's effect on mitochondria. This means we are looking for an experiment that is either irrelevant to mitochondrial function or is based on a flawed biological premise.

### Step 2: Analyze Each Experimental Option

**A) Confocal fluorescence microscopy after Mito-RTP staining of the cells**
*   **What it measures:** This experiment uses a fluorescent probe named "Mito-RTP". The "Mito-" prefix indicates that the probe is designed to specifically accumulate within the mitochondria. Probes like this are used to measure key parameters of mitochondrial function or status, such as temperature (as a proxy for metabolic rate) or the presence of reactive oxygen species (ROS).
*   **Relevance:** Changes in mitochondrial metabolism, a core function, directly lead to changes in heat production and ROS levels. Visualizing these changes with high-resolution microscopy is a valid and direct method to observe if the drug affects mitochondrial activity or induces stress.
*   **Careful Point:** This is a direct and specific method for investigating a drug's effect on the internal state of mitochondria. This experiment is helpful.

**B) Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant**
*   **What it measures:** This describes a luciferase-based assay, which measures the concentration of ATP. The enzyme luciferase uses ATP to produce light from its substrate, luciferin.
*   **Relevance:** Mitochondria are the primary producers of ATP in the cell. Measuring ATP levels is a fundamental way to assess the primary output of mitochondrial function. A drug that impairs mitochondria would almost certainly lead to a change in cellular ATP levels.
*   **Careful Point:** The prompt specifies measuring ATP in the *supernatant* (extracellular ATP). While measuring intracellular ATP is a more direct measure of production, a significant release of ATP into the supernatant is a known indicator of cell stress or death, which can be triggered by mitochondrial toxins. Therefore, this is a valid, albeit indirect, method to screen for drug-induced effects related to cellular energy and viability. This experiment is helpful.

**C) Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide**
*   **What it measures:** The long chemical name refers to the common fluorescent dye **JC-1**. JC-1 is used to measure the mitochondrial membrane potential (ΔΨm), a critical indicator of mitochondrial health.
*   **Relevance:** The mitochondrial membrane potential is directly generated by the redox reactions of the electron transport chain and is necessary for ATP synthesis. A high potential indicates healthy, functioning mitochondria. A drop in potential is a classic sign of mitochondrial dysfunction. Flow cytometry allows for a rapid, quantitative measurement of this change across a large population of cells.
*   **Careful Point:** This is a standard, powerful, and direct method for assessing a drug's impact on mitochondrial function. This experiment is helpful.

**D) Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit**
*   **What it measures:** This experiment first isolates mitochondria from the cell. Then, it attempts to measure the rate of glucose uptake by these isolated mitochondria.
*   **Relevance:** This is where the critical flaw lies. The process of glucose uptake from the environment is carried out by glucose transporters (like GLUTs) located on the **cell's plasma membrane**, not on the mitochondrial membranes. The initial breakdown of glucose (glycolysis) occurs in the **cytoplasm**. Mitochondria import the products of glycolysis (primarily pyruvate), not glucose itself.
*   **Careful Point:** The experimental design is based on a false biological premise. Applying a glucose uptake assay to an isolated mitochondrial fraction is a nonsensical procedure because the organelle lacks the necessary machinery to perform the function being measured. The experiment will not yield any meaningful data about mitochondrial function. This experiment is **not helpful**.

### Final Decision
Experiments A, B, and C all describe valid methods to measure key parameters directly or indirectly related to mitochondrial function: metabolic activity/stress, ATP levels/cell viability, and membrane potential, respectively. Experiment D is based on a fundamental biological error and is therefore the only experiment that will not help discover the mitochondrial role of the drug.

<<<D>>>
"""

# Run the check
result = check_answer(question, llm_answer)
print(result)