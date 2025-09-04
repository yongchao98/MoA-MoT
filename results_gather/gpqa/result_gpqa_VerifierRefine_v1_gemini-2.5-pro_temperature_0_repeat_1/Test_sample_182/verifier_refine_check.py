import re

def calculate_ihd(formula: str) -> float:
    """
    Calculates the Index of Hydrogen Deficiency (IHD) for a given molecular formula.
    The formula for IHD is: C - H/2 - X/2 + N/2 + 1
    """
    # Find all element-count pairs in the formula string
    atoms = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    
    counts = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'X': 0}
    
    for element, count_str in atoms:
        count = int(count_str) if count_str else 1
        if element == 'C':
            counts['C'] += count
        elif element == 'H':
            counts['H'] += count
        elif element == 'N':
            counts['N'] += count
        elif element in ['F', 'Cl', 'Br', 'I']:
            counts['X'] += count
        # Oxygen does not affect the IHD calculation
        
    C = counts['C']
    H = counts['H']
    N = counts['N']
    X = counts['X']
    
    ihd = C - (H / 2) - (X / 2) + (N / 2) + 1
    return ihd

def check_answer():
    """
    Checks the correctness of the provided answer for the chemistry question.
    """
    # --- Step 1: Define reactant and product based on the problem statement ---
    
    # Reactant: 2-formyl-5-vinylcyclohex-3-enecarboxylic acid
    # Derived formula: C10H12O3
    reactant_formula = "C10H12O3"
    
    # Reaction: with red phosphorus and excess of HI.
    # This reaction reduces all C=C and C=O bonds and removes all oxygen atoms,
    # resulting in a saturated monocyclic hydrocarbon.
    # The number of carbon atoms remains the same (10).
    # The general formula for a saturated monocyclic hydrocarbon is CnH2n.
    # Therefore, the product formula is C10H20.
    product_formula = "C10H20"
    
    # The answer provided by the LLM
    llm_answer = 1

    # --- Step 2: Calculate IHD for verification ---
    
    # IHD of the reactant (for completeness)
    reactant_ihd = calculate_ihd(reactant_formula)
    # Expected reactant IHD: 1 ring + 1(C=C in ring) + 1(C=C in vinyl) + 1(C=O in formyl) + 1(C=O in acid) = 5
    if reactant_ihd != 5:
        return f"Error in reactant analysis. Calculated reactant IHD is {reactant_ihd}, but it should be 5."

    # IHD of the product
    product_ihd = calculate_ihd(product_formula)
    
    # --- Step 3: Compare calculated IHD with the LLM's answer ---
    
    print(f"Reactant Formula: {reactant_formula}")
    print(f"Calculated IHD of Reactant: {reactant_ihd}")
    print("-" * 30)
    print(f"Product Formula: {product_formula}")
    print(f"Calculated IHD of Product: {product_ihd}")
    print(f"LLM's Answer for Product IHD: {llm_answer}")
    print("-" * 30)

    if product_ihd == llm_answer:
        # The product is a saturated monocyclic alkane (a substituted cyclohexane).
        # It has 1 ring and 0 pi bonds.
        # Therefore, its IHD must be 1.
        # The calculation confirms this.
        return "Correct"
    else:
        return (f"Incorrect. The calculated IHD of the product ({product_formula}) is {product_ihd}, "
                f"but the provided answer is {llm_answer}. The product is a saturated monocyclic alkane, "
                f"which must have an IHD of 1 (for the one ring).")

# Run the check
result = check_answer()
print(result)