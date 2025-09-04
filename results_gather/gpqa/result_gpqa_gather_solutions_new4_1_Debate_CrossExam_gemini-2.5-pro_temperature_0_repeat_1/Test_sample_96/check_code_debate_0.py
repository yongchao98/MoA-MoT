import re

def check_correctness(question, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the given biological question.

    Args:
        question (str): The question text.
        llm_answer_text (str): The full text of the LLM's answer, including the final choice.

    Returns:
        str: "Correct" if the answer is correct, otherwise a reason for the error.
    """

    # Extract the final answer choice from the format <<<X>>>
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: The answer format is invalid. It should end with '<<<X>>>' where X is A, B, C, or D."
    
    llm_choice = match.group(1)

    # Define the correct answer and the scientific reasoning
    correct_choice = 'D'
    
    scientific_facts = {
        'A': "Incorrect. Chiasmata resolution is a meiotic process. Errors here can cause aneuploidy but do not explain the difference in phenotypic severity after fertilization.",
        'B': "Incorrect. Spindle attachment is a meiotic process. Errors here can cause aneuploidy but do not explain the difference in phenotypic severity after fertilization.",
        'C': "Incorrect. DNA polymerase alpha is involved in general DNA replication, not the specific dosage compensation mechanism that mitigates Klinefelter's syndrome.",
        'D': "Correct. The less severe phenotype in Klinefelter's syndrome is due to X-chromosome inactivation (XCI). This is a post-zygotic epigenetic process where one X chromosome is silenced. A key molecular event in this process is chromatin methylation by histone methyltransferases."
    }

    # Check if the LLM's choice matches the correct choice
    if llm_choice == correct_choice:
        # Further check if the reasoning provided in the LLM's text is sound.
        # The provided reasoning correctly identifies XCI, epigenetic modification, and the post-zygotic timing.
        reasoning_text = llm_answer_text.lower()
        required_keywords = [
            'x-chromosome inactivation', 'dosage compensation', 'epigenetic', 
            'barr body', 'histone methyltransferase', 'post-zygote'
        ]
        
        # Check if the core concepts are present in the explanation
        if all(keyword in reasoning_text for keyword in required_keywords):
            return "Correct"
        else:
            # The choice is correct, but the reasoning might be incomplete.
            # For this problem, the choice itself is the primary check.
            return "Correct"
    else:
        # The LLM's choice is incorrect. Provide a detailed reason.
        reason_for_error = f"Incorrect. The selected answer was {llm_choice}, but the correct answer is {correct_choice}.\n"
        reason_for_error += f"Reason for {llm_choice}: {scientific_facts[llm_choice]}\n"
        reason_for_error += f"Reason for {correct_choice}: {scientific_facts[correct_choice]}"
        return reason_for_error

# --- Execution of the check ---
question_text = """
In gametogenesis, the germline undergoes two cycles of meiosis to give rise to haploid gametes. The goal of meiosis, unlike mitosis, is to separate not only the replicated sister chromatids but also the homologous chromosome pairs to form haploid cells. Fertilized gametes enter embryogenesis and ultimately develop into an organism. Meiosis is not always perfectly executed, and most such cases result in embryonic lethality. Viable events can already be detected in the karyotype of embryonic cells from the amniotic fluid. One such event results in Klinefelter's syndrome in boys. Which molecular mechanism is responsible for less prominent phenotypic consequences of Klinefelter's syndrome compared to better known Down's syndrome?

A) chiasmata resolution by separase in diakinesis
B) attachment of spindle to kinetochores in the metaphase I
C) progression of the polymerase alpha in the morula/blastocyst
D) chromatin methylation by histone methyltransferases in the post-zygote
"""

# This is the final answer provided by the user's prompt
final_llm_answer = """
Here is a step-by-step analysis of the question based on the provided candidate answers.

### Step 1: Deconstruct the Question
The core of the question is to identify the molecular reason why having an extra X chromosome (Klinefelter's syndrome, 47,XXY) results in a less severe phenotype than having an extra chromosome 21 (Down's syndrome, Trisomy 21). This requires comparing the biological response to a sex chromosome aneuploidy versus an autosomal aneuploidy.

### Step 2: Analyze the Consensus from Candidate Answers
A vast majority of the candidate answers (at least 15 out of 17) converge on a single explanation. The key points from this consensus are:
*   **Dosage Compensation:** The fundamental reason for the difference in severity is a mechanism called **X-chromosome inactivation (XCI)**, also known as Lyonization.
*   **Mechanism of XCI:** In any mammalian cell with more than one X chromosome, all but one are epigenetically silenced to ensure the "dose" of X-linked genes is appropriate. This process occurs in typical XX females and also in individuals with Klinefelter's syndrome (XXY).
*   **Consequence for Klinefelter's:** In an XXY individual, one X chromosome is inactivated and becomes a condensed, mostly silent structure called a **Barr body**. This largely normalizes the gene expression from the X chromosomes, thus mitigating the potential harm of having an extra one.
*   **Contrast with Down's Syndrome:** There is no equivalent, large-scale inactivation mechanism for autosomes. Therefore, in an individual with Down's syndrome, all three copies of chromosome 21 remain active, leading to a significant (~50%) overexpression of hundreds of genes. This massive gene dosage imbalance causes the more severe and wide-ranging phenotype.

### Step 3: Evaluate the Multiple-Choice Options Based on the Consensus
The question asks for the *molecular mechanism* responsible for this difference. Based on the analysis above, we are looking for an option that describes a key molecular event within X-chromosome inactivation.

*   **A) chiasmata resolution by separase in diakinesis** and **B) attachment of spindle to kinetochores in the metaphase I:** These options describe processes that occur during meiosis. Errors in these processes (nondisjunction) are the *cause* of aneuploidies like Klinefelter's and Down's syndromes. They do not, however, explain the difference in the *consequences* or severity of these syndromes after fertilization. All candidate answers correctly dismiss these options.

*   **C) progression of the polymerase alpha in the morula/blastocyst:** DNA polymerase alpha is an enzyme essential for DNA replication, a general process required for all cell division in the early embryo (morula/blastocyst stage). This is not a specific mechanism for gene silencing or dosage compensation that would differentiate the outcomes of these two syndromes.

*   **D) chromatin methylation by histone methyltransferases in the post-zygote:** This option describes a specific molecular event. X-chromosome inactivation is an epigenetic process that occurs after fertilization ("post-zygote"). A central part of this process is the modification of chromatin to make it compact and transcriptionally silent. This involves recruiting enzymes, including **histone methyltransferases**, which add methyl groups to histone proteins (a form of **chromatin methylation**). This is a hallmark of silenced chromatin and is essential for establishing and maintaining the inactive state of the Barr body. This option accurately describes a core component of the mechanism that lessens the severity of Klinefelter's syndrome.

### Step 4: Final Conclusion
The less prominent phenotype of Klinefelter's syndrome is a direct result of X-chromosome inactivation, which compensates for the extra X chromosome. The molecular basis of this inactivation is epigenetic silencing. Option D provides the most accurate description of a key molecular event in this process among the choices provided. The overwhelming majority of the candidate answers, through detailed and consistent reasoning, support this conclusion.

<<<D>>>
"""

result = check_correctness(question_text, final_llm_answer)
print(result)