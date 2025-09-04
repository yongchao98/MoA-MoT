import re

def check_answer(question_text, llm_answer_text):
    """
    Checks the correctness of an LLM's answer to a molecular biology question.

    The function simulates the logical deduction process based on the facts
    provided in the question. It verifies if the chosen answer is the most
    plausible explanation for the observed experimental outcome.
    """

    # --- Step 1: Extract facts and constraints from the question ---

    # Fact 1: The core molecular mechanism involves Cre-Lox recombination.
    # The key element is the lox2272 site left after recombination.
    # This is a known biological constant.
    lox_site_length = 34  # in base pairs

    # Fact 2: The genetic code is read in triplets (codons).
    codon_length = 3

    # Fact 3: The construct uses a strong, ubiquitous promoter.
    promoter_is_strong = "CBA promoter" in question_text

    # Fact 4: An in-vitro test confirmed protein expression before recombination.
    in_vitro_test_successful = "test the construct expression on the Western blot after plasmid transfection" in question_text

    # Fact 5: The observed outcome in the mouse model.
    observed_outcome = "You do not observe a green signal"
    
    # --- Step 2: Define the options and the LLM's answer ---

    # The provided analysis block establishes a canonical mapping for the options,
    # even if the candidate answers have inconsistent lettering.
    options = {
        'A': "ligand and the receptor are in a paracrine relationship",
        'B': "the receptor and the eGFP are not in the frame",
        'C': "the enhancer for the ligand and receptor expression is missing",
        'D': "the receptor-eGFP construct is stuck in the Golgi"
    }

    # Extract the letter from the LLM's final answer format, e.g., <<<B>>>
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Invalid answer format. Expected '<<<X>>>' where X is A, B, C, or D."
    
    llm_choice = match.group(1)
    
    # --- Step 3: Analyze each option based on the facts ---

    analysis = {}
    
    # Analyze Option A: Paracrine relationship
    # This describes a biological function between cells, not a molecular reason
    # for protein synthesis failure within a cell.
    analysis['A'] = {
        'is_plausible': False,
        'reason': "This describes a biological function, which is irrelevant to the technical failure of protein synthesis from the reporter construct."
    }

    # Analyze Option B: Frameshift
    # This is the core technical check.
    is_frameshift = (lox_site_length % codon_length) != 0
    if is_frameshift:
        analysis['B'] = {
            'is_plausible': True,
            'reason': f"The residual {lox_site_length} bp lox2272 site is not divisible by the codon length of {codon_length}, causing a frameshift mutation. This prevents the synthesis of a functional eGFP protein, perfectly explaining the absence of a green signal."
        }
    else:
        analysis['B'] = {
            'is_plausible': False,
            'reason': f"A {lox_site_length} bp lox site would not cause a frameshift, so this cannot be the reason."
        }

    # Analyze Option C: Missing enhancer
    # This is contradicted by the problem statement.
    if promoter_is_strong and in_vitro_test_successful:
        analysis['C'] = {
            'is_plausible': False,
            'reason': "This is contradicted by the use of the strong, ubiquitous CBA promoter and the successful Western blot test, which confirmed the construct's ability to express protein."
        }
    else:
        analysis['C'] = {
            'is_plausible': True, # Should not happen with the given text
            'reason': "This might be plausible if the promoter was weak or untested."
        }

    # Analyze Option D: Stuck in Golgi
    # This contradicts the observed outcome.
    if "no green signal" in observed_outcome:
        analysis['D'] = {
            'is_plausible': False,
            'reason': "A protein stuck in the Golgi would still be fluorescent, leading to a mislocalized signal. The observation was a complete *absence* of signal, which points to a synthesis failure, not a trafficking issue."
        }
    else:
        analysis['D'] = {
            'is_plausible': True, # Should not happen with the given text
            'reason': "This might be plausible if the observation was a mislocalized signal."
        }

    # --- Step 4: Determine the correct answer and check against the LLM's choice ---
    
    correct_option = None
    for option, result in analysis.items():
        if result['is_plausible']:
            correct_option = option
            break
            
    if llm_choice == correct_option:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer '{llm_choice}' is wrong.\n"
                f"Reason why '{llm_choice}' is wrong: {analysis[llm_choice]['reason']}\n"
                f"The correct answer is '{correct_option}'.\n"
                f"Reason why '{correct_option}' is correct: {analysis[correct_option]['reason']}")


# The question text and the final answer provided by the user
question = """
You identified a new ligand-receptor pair initially expressed in the neural crest cells and want to understand the developmental fate of their interaction. The receptor is a monomeric transmembrane barrel with eight alpha helices, four extracellular helices, and three beta sheets that create a binding site for the ligand. The ligand is a coiled-coil of two alpha helices. You create a bicistronic lox-Cre vector with the ligand ORF at the 5' position, fused to mCherry; the receptor ORF is at the 3' and fused to the eGFP at the C-terminus. The ligand is under the CBA promoter, and upstream of the receptor is the IRES regulatory element. Between the ligand ORF and mCherry, there is a stop in a loxP-stop-loxP cassette, and between the receptor and the eGFP, there is a stop in the lox2272-stop-lox2272 cassette. You had also raised antibodies against the overexpressed proteins and now use them to test the construct expression on the Western blot after plasmid transfection to primary astrocyte culture. Total proteins give a signal proportional to the actin loading. You finally breed a mouse model homozygous for your construct in the wt background. You cross it to SOX10-Cre hemizygous mice; after the offspring is born you visualize them under the fluorescent confocal microscope. You do not observe a green signal. Which of the answers below explains the most likely reason for it?
"""

llm_answer = """
My plan is to first confirm the length of the lox2272 site, as this is the most critical piece of information for solving the problem. Then, I will use this information to analyze the molecular consequences of the Cre-recombination event and evaluate each of the possible answers.

I will now execute the first step of my plan.
<search>lox2272 site length</search>
Here are the search results.
<information>
[1] The loxP site is a 34 bp sequence consisting of two 13 bp palindromic repeats separated by an 8 bp spacer region. The lox2272 site is a variant of the loxP site with mutations in the spacer region. It is still 34 bp long.

[2] Cre-Lox recombination is a powerful tool for genetic manipulation. When a DNA segment flanked by two lox sites is excised, a single 34 bp lox site remains as a "scar". If this scar is within a coding sequence, it can cause a frameshift mutation unless the reading frame is specifically corrected, because 34 is not divisible by 3.

[3] To create in-frame fusions using Cre-Lox, engineers must use modified lox sites or add/remove bases adjacent to the lox site to ensure the total number of base pairs in the scar is a multiple of three. Standard loxP, lox2272, and loxN sites are all 34 bp and will cause a +1 frameshift if placed directly between two coding sequences.
</information>
The search confirms that the `lox2272` site is 34 base pairs long. This is the critical piece of information needed to solve the problem. I can now proceed with the analysis.

**Step-by-step Analysis**

1.  **Deconstruct the Experimental Design:** The experiment uses a Cre-inducible system to create a Receptor-eGFP fusion protein. The relevant part of the construct is `Receptor ORF -> lox2272-stop-lox2272 -> eGFP ORF`.
2.  **Analyze the Molecular Event:** In the `SOX10-Cre` mouse, Cre recombinase is present in neural crest cells. Its function is to recognize the two `lox2272` sites and excise the intervening "stop" cassette.
3.  **Identify the Consequence of Recombination:** This recombination event leaves a single `lox2272` site as a "scar" in the DNA, positioned between the Receptor's coding sequence and the eGFP's coding sequence.
4.  **Identify the Critical Flaw:** As confirmed by the search, the `lox2272` scar is **34 base pairs (bp)** long. The genetic code is read in triplets (codons). For a fusion protein to be synthesized correctly, the number of base pairs in any sequence linking the two coding regions must be a multiple of 3.
5.  **Determine the Result:** Since 34 is not divisible by 3 (34 ÷ 3 = 11 with a remainder of 1), the 34 bp scar introduces a **frameshift mutation**. The ribosome translates the receptor protein correctly, but its reading frame is shifted when it encounters the scar. It then begins to read the eGFP sequence in the wrong frame, producing a nonsensical amino acid sequence and almost certainly encountering a premature stop codon. No functional, fluorescent eGFP protein can be synthesized. This perfectly explains the complete absence of a green signal.

**Evaluation of the Answer Choices**

*   **A) ligand and the receptor are in a paracrine relationship:** This describes a biological function (how cells signal to each other) and is irrelevant to the technical failure of protein synthesis from the reporter construct within a single cell.
*   **B) the receptor and the eGFP are not in the frame:** This is the correct explanation. The frameshift mutation caused by the 34 bp lox2272 scar prevents the correct translation of the eGFP protein.
*   **C) the enhancer for the ligand and receptor expression is missing:** This is incorrect. The construct uses the strong, ubiquitous CBA promoter, and the initial Western blot in astrocytes confirmed that the promoter and IRES are functional, allowing the non-fused receptor protein to be expressed.
*   **D) the receptor-eGFP construct is stuck in the Golgi:** This is unlikely. If the protein were synthesized but mis-trafficked, it would still be fluorescent. A confocal microscope would detect a green signal, just mislocalized to the Golgi apparatus. The observation is a complete *absence* of a signal, which points to a failure in synthesis, not trafficking.

Based on this analysis, the frameshift mutation is the only explanation that accounts for the complete and specific failure to produce a green signal.

<<<B>>>
"""

# Run the check
result = check_answer(question, llm_answer)
print(result)