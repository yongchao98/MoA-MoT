import re

def check_correctness_of_answer(question: str, final_answer: str) -> str:
    """
    Checks the correctness of the final answer for the given biological question.

    The function programmatically applies constraints derived from the question to determine the correct option
    and then compares it with the provided final answer.

    Args:
        question: The full text of the question, including the options.
        final_answer: The final answer provided by the LLM, including the reasoning and the <<<X>>> tag.

    Returns:
        A string indicating "Correct" or a detailed reason for the incorrectness.
    """
    # Extract the letter from the final answer format, e.g., <<<D>>>
    answer_match = re.search(r'<<<([A-D])>>>', final_answer)
    if not answer_match:
        return "Invalid answer format: The final answer must end with '<<<X>>>' where X is A, B, C, or D."
    
    provided_answer_letter = answer_match.group(1)

    # Define the options and their biological properties based on the question's context.
    # This allows for programmatic evaluation.
    options = {
        'A': {
            'text': 'chiasmata resolution by separase in diakinesis',
            'timing': 'meiotic',  # Occurs during meiosis, before fertilization.
            'function': 'cause_of_aneuploidy' # An error here causes the syndrome.
        },
        'B': {
            'text': 'progression of the polymerase alpha in the morula/blastocyst',
            'timing': 'post-zygotic', # Occurs after fertilization.
            'function': 'general_cell_process' # DNA replication is not specific to dosage compensation.
        },
        'C': {
            'text': 'attachment of spindle to kinetochores in the metaphase I',
            'timing': 'meiotic', # Occurs during meiosis, before fertilization.
            'function': 'cause_of_aneuploidy' # An error here causes the syndrome.
        },
        'D': {
            'text': 'chromatin methylation by histone methyltransferases in the post-zygote',
            'timing': 'post-zygotic', # Occurs after fertilization.
            'function': 'mitigation_of_aneuploidy' # This is the specific mechanism (X-inactivation) that lessens the severity.
        }
    }

    # --- Constraint 1: The mechanism must explain the difference in phenotypic severity, not the cause of the syndrome. ---
    # The question asks why Klinefelter's is less severe, which implies a post-fertilization (post-zygotic) mechanism
    # that mitigates the effects. Mechanisms that cause the syndrome (meiotic errors) are incorrect.
    
    # Filter for options that describe a mitigating function.
    valid_options = {key: val for key, val in options.items() if val['function'] == 'mitigation_of_aneuploidy'}
    
    # If we apply the timing constraint as well, the result is the same.
    # valid_options = {key: val for key, val in options.items() if val['timing'] == 'post-zygotic' and val['function'] != 'general_cell_process'}

    if len(valid_options) != 1:
        # This case should not be reached with the current options but is a safeguard.
        return f"Error in checking logic: Expected 1 correct option based on constraints, but found {len(valid_options)}."

    correct_answer_letter = list(valid_options.keys())[0]

    # --- Final Verification ---
    if provided_answer_letter == correct_answer_letter:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The provided answer is <<<{provided_answer_letter}>>>.\n\n"
            f"The question asks for the mechanism that makes Klinefelter's syndrome (extra X chromosome) *less severe* than Down's syndrome (extra autosome). This requires identifying a process that *mitigates* the effect of the extra chromosome, not one that *causes* it.\n\n"
            f"Constraint 1: The mechanism must occur post-fertilization (post-zygotic) to affect the developing phenotype. Options A and C describe meiotic events, which are potential *causes* of the syndrome, not mitigators of its effects. They are therefore incorrect.\n\n"
            f"Constraint 2: The mechanism must be specific to the differential handling of a sex chromosome versus an autosome. Option B (DNA replication) is a general cellular process and does not explain the specific dosage compensation for the X chromosome.\n\n"
            f"Constraint 3: The correct mechanism is X-chromosome inactivation (XCI), an epigenetic process that silences the extra X chromosome. Option D, '{options['D']['text']}', accurately describes a key molecular event of XCI.\n\n"
            f"Therefore, the only option that satisfies all constraints is D. The provided answer <<<{provided_answer_letter}>>> is wrong."
        )
        return reason

# The user's provided question and final answer to be checked.
question_text = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
In gametogenesis, the germline undergoes two cycles of meiosis to give rise to haploid gametes. The goal of meiosis, unlike mitosis, is to separate not only the replicated sister chromatids but also the homologous chromosome pairs to form haploid cells. Fertilized gametes enter embryogenesis and ultimately develop into an organism. Meiosis is not always perfectly executed, and most such cases result in embryonic lethality. Viable events can already be detected in the karyotype of embryonic cells from the amniotic fluid. One such event results in Klinefelter's syndrome in boys. Which molecular mechanism is responsible for less prominent phenotypic consequences of Klinefelter's syndrome compared to better known Down's syndrome?

A) chiasmata resolution by separase in diakinesis
B) progression of the polymerase alpha in the morula/blastocyst
C) attachment of spindle to kinetochores in the metaphase I
D) chromatin methylation by histone methyltransferases in the post-zygote

Once you feel you are ready for the final answer, directly return the answer with the format <<<answer content>>> at the end of your response, e.g. <<<C>>>, <<<A>>>
"""

final_answer_from_llm = """
Here is a step-by-step analysis to determine the final answer based on the provided candidate responses.

### Step 1: Deconstruct the Core Question
The question asks for the specific molecular mechanism that explains why Klinefelter's syndrome (47,XXY) has less severe phenotypic consequences than Down's syndrome (Trisomy 21). The key is to identify a process that *mitigates* the effects of the extra chromosome, not the process that *causes* the aneuploidy in the first place.

### Step 2: Identify the Key Biological Principle
A thorough analysis of the candidate answers reveals a strong and accurate consensus on the underlying biological principle:
*   **The Difference:** The core difference is that Klinefelter's involves an extra **sex chromosome** (X), while Down's involves an extra **autosome** (chromosome 21).
*   **The Mechanism:** The body has a natural process for handling extra X chromosomes called **X-chromosome inactivation (XCI)**, or Lyonization. In this process, all but one X chromosome in a cell are epigenetically silenced.
*   **Application to Syndromes:**
    *   In Klinefelter's (XXY), one X chromosome is inactivated, leaving one active X and one Y. This largely normalizes the gene dosage, making the phenotype much less severe.
    *   In Down's syndrome, there is no equivalent mechanism to inactivate an autosome. All three copies of chromosome 21 remain active, leading to a significant gene dosage imbalance and a more severe phenotype.

### Step 3: Analyze the Molecular Basis of X-Chromosome Inactivation
The candidate answers correctly describe XCI as an **epigenetic** process that occurs **post-zygotically** (after fertilization, in the early embryo). A central molecular event in this process is the modification of chromatin (the complex of DNA and histone proteins) to make it compact and transcriptionally silent. This involves:
*   **Chromatin methylation by histone methyltransferases**: These enzymes add methyl groups to histone proteins, which is a key repressive mark that signals for gene silencing.

### Step 4: Evaluate the Canonical Options
Let's evaluate the options provided in the original prompt against this biological understanding.

*   **A) chiasmata resolution by separase in diakinesis:** This is a process during meiosis. An error here can *cause* aneuploidy but does not explain the difference in severity *after* fertilization. This is a distractor related to the cause, not the consequence.
*   **B) progression of the polymerase alpha in the morula/blastocyst:** This refers to DNA replication, a general process required for all cell division. It is not the specific mechanism that silences an entire chromosome for dosage compensation.
*   **C) attachment of spindle to kinetochores in the metaphase I:** Similar to option A, this is a meiotic process where errors can *cause* aneuploidy. It does not explain the downstream phenotypic differences.
*   **D) chromatin methylation by histone methyltransferases in the post-zygote:** This option perfectly describes a core molecular event of X-chromosome inactivation. It correctly identifies the epigenetic mechanism ("chromatin methylation"), the enzymes involved ("histone methyltransferases"), and the timing ("in the post-zygote"). This is the direct explanation for the mitigated phenotype in Klinefelter's syndrome.

### Step 5: Synthesize and Conclude
The reasoning across almost all candidate answers is sound and points directly to X-chromosome inactivation. The wide variation in the final letter chosen (A, B, C, D) strongly suggests that the options were shuffled for each agent. By ignoring the agents' final letter choices and focusing on their high-quality, consensus reasoning, we can map that reasoning to the fixed set of options in the prompt. The consensus reasoning aligns perfectly and exclusively with option D.

<<<D>>>
"""

# Execute the check
result = check_correctness_of_answer(question_text, final_answer_from_llm)
print(result)