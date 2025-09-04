import re

class ChemistryChecker:
    """
    This class checks the correctness of a multi-step chemistry problem solution.
    It verifies the reaction pathway and the symmetry analysis of the final product.
    """

    def __init__(self, question, llm_answer_text):
        """
        Initializes the checker with the question and the LLM's answer.
        """
        self.question = question
        self.llm_answer_text = llm_answer_text
        # This dictionary stores the correct point groups for potential products.
        # This is a standard way to verify chemical reasoning without running complex simulations.
        self.point_groups = {
            "p-nitrobenzoate anion": "C2v",
            "(E)-4-(4-nitrophenyl)but-3-en-2-one": "Cs",
            "(E)-4,4'-azodibenzoic acid": "C2h",
            "(E)-4,4'-azoxybis(benzoic acid)": "Cs",
            "1,5-bis(4-nitrophenyl)penta-1,4-dien-3-one": "C2v"
        }
        self.options_map = {'A': 'c3', 'B': 'd2h', 'C': 'c2h', 'D': 'cs'}

    def check_correctness(self):
        """
        Performs a step-by-step check of the LLM's reasoning and final answer.
        """
        # 1. Check the final answer format and value
        match = re.search(r'<<<([A-D])>>>', self.llm_answer_text)
        if not match:
            return "Format Error: The final answer is not in the required '<<<A>>>' format."
        
        final_choice = match.group(1)
        if final_choice != 'D':
            return f"Incorrect Answer: The final choice is '{final_choice}', but the correct choice is 'D'."

        # 2. Verify the reaction pathway identified in the reasoning
        # Step 1: Nitration of toluene
        if "p-nitrotoluene" not in self.llm_answer_text.lower():
            return "Reasoning Error: Step 1 is incorrect. The major product of toluene nitration is p-nitrotoluene."
        
        # Step 2: Oxidation of p-nitrotoluene
        # The most logical intermediate for the next step is the aldehyde.
        if "p-nitrobenzaldehyde" not in self.llm_answer_text.lower():
            return "Reasoning Error: Step 2 is incorrect. The most logical intermediate for the subsequent Claisen-Schmidt condensation is p-nitrobenzaldehyde."

        # Step 3: Claisen-Schmidt condensation
        product3_name = "(E)-4-(4-nitrophenyl)but-3-en-2-one"
        # Normalize the name for checking
        normalized_product3_name = product3_name.lower().replace("(e)-", "").replace("-", "")
        normalized_llm_text = self.llm_answer_text.lower().replace("(e)-", "").replace("-", "")
        if normalized_product3_name not in normalized_llm_text:
            return f"Reasoning Error: Step 3 is incorrect. The product of the Claisen-Schmidt condensation should be {product3_name}."

        # 3. Verify the symmetry analysis
        correct_point_group = self.point_groups.get(product3_name)
        if not correct_point_group:
             return f"Internal Check Error: Point group for '{product3_name}' not defined."

        # Check if the LLM correctly identifies the point group in its text
        if f"point group is {correct_point_group.lower()}" not in self.llm_answer_text.lower() and \
           f"symmetry is {correct_point_group.lower()}" not in self.llm_answer_text.lower():
            return f"Symmetry Error: The reasoning fails to correctly identify the point group of {product3_name} as {correct_point_group}."

        # 4. Verify the final choice matches the reasoning
        if self.options_map.get(final_choice).lower() != correct_point_group.lower():
            return f"Answer Mismatch: The reasoning correctly identifies the point group as {correct_point_group}, but the selected option '{final_choice}' is incorrect."

        return "Correct"

# The user's provided question and the LLM's answer to be checked.
question = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
toluene is treated with nitric acid and sulfuric acid, forming product 1.

1 is treated with MnO2 and H2SO4, forming product 2.

2 is treated with acetone and aqueous sodium hydroxide, forming product 3.

what is the molecular symmetry group of 3?

A) c3
B) d2h
C) c2h
D) cs
"""

llm_answer = """
Here is a step-by-step analysis of the problem, considering the various interpretations presented by the candidate answers.

### **Step 1: Toluene to Product 1**

*   **Reaction:** Toluene is treated with nitric acid (HNO₃) and sulfuric acid (H₂SO₄).
*   **Analysis:** This is a standard electrophilic aromatic substitution, specifically nitration. The methyl group (-CH₃) is an ortho-, para-director. Due to steric hindrance, the major product is the para isomer.
*   **Conclusion:** All candidate answers correctly identify **Product 1 as p-nitrotoluene**.

### **Step 2: Product 1 to Product 2**

*   **Reaction:** p-nitrotoluene is treated with manganese dioxide (MnO₂) and sulfuric acid (H₂SO₄).
*   **Analysis:** This is an oxidation of the benzylic methyl group. There are two main possibilities debated among the candidates:
    1.  **Oxidation to Aldehyde:** Formation of p-nitrobenzaldehyde.
    2.  **Oxidation to Carboxylic Acid:** Formation of p-nitrobenzoic acid.
*   **Careful Point:** The choice between these two products is critical as it dictates the next step. In multi-step synthesis problems, the identity of an intermediate is often clarified by the subsequent reaction. The reagents in Step 3 (acetone and NaOH) are the classic setup for a Claisen-Schmidt condensation, which requires an aldehyde (or ketone) as a reactant. A carboxylic acid would simply undergo an acid-base reaction. Therefore, the most logical intended product is the aldehyde.
*   **Conclusion:** **Product 2 is p-nitrobenzaldehyde**.

### **Step 3: Product 2 to Product 3**

*   **Reaction:** p-nitrobenzaldehyde is treated with acetone and aqueous sodium hydroxide (NaOH).
*   **Analysis:** This is a base-catalyzed crossed-aldol condensation (specifically, a Claisen-Schmidt condensation). The base (NaOH) deprotonates acetone to form a nucleophilic enolate. The enolate attacks the carbonyl of p-nitrobenzaldehyde. The resulting intermediate rapidly dehydrates to form a stable, conjugated product.
*   **Careful Point:** Acetone has α-protons on both sides of the carbonyl. While a double condensation is possible, the simpler and most common product shown in textbook examples is the single condensation product.
*   **Conclusion:** **Product 3 is (E)-4-(4-nitrophenyl)but-3-en-2-one**.

### **Step 4: Molecular Symmetry of Product 3**

*   **Structure:** (E)-4-(4-nitrophenyl)but-3-en-2-one.
    `O₂N-C₆H₄-CH=CH-C(=O)-CH₃`
*   **Analysis:** To determine the point group, we must identify the symmetry elements.
    *   The molecule has an extended conjugated system, causing the heavy-atom framework (C, N, O) to be essentially planar in its most stable conformation.
    *   **Plane of Symmetry (σ):** The plane containing this planar framework is a mirror plane. The molecule is reflected onto itself by this plane.
    *   **Axis of Rotation (Cₙ, n>1):** There is no C₂ axis or any higher-order rotational axis. The two ends of the molecule are different (a p-nitrophenyl group vs. an acetyl group), so no rotation can leave the molecule unchanged.
    *   **Center of Inversion (i):** There is no center of inversion.
*   **Conclusion:** The molecule possesses only the identity element (E) and a single plane of symmetry (σ). A molecule with only these elements belongs to the **Cₛ** point group.

### **Final Decision**

The most chemically sound and logical reaction sequence, designed to test standard organic reactions, leads to (E)-4-(4-nitrophenyl)but-3-en-2-one. This molecule has **Cₛ** symmetry. According to the options provided in the prompt:
A) c3
B) d2h
C) c2h
D) cs

The correct option is D.

<<<D>>>
"""

# Instantiate the checker and run the check
checker = ChemistryChecker(question, llm_answer)
result = checker.check_correctness()
print(result)