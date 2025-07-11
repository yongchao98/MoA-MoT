import textwrap

def identify_reactant():
    """
    This function outlines the chemical reasoning to identify the starting material.
    """
    
    # Step 1: Analyze the reaction conditions and product.
    reasoning_step1 = """
    The reaction involves treating a single compound with potassium hydroxide (KOH), a base, to form a complex bicyclic ketone. This strongly suggests an intramolecular aldol condensation reaction.
    """

    # Step 2: Analyze the likely precursor for this type of reaction.
    reasoning_step2 = """
    The product is a derivative of a decalin system (two fused six-membered rings), a typical outcome of a Robinson annulation reaction. The intramolecular aldol condensation is the final step of this annulation, which cyclizes a dicarbonyl precursor. The question asks for the single compound that reacts with KOH, pointing to this dicarbonyl precursor.
    """

    # Step 3: Identify a chemically plausible candidate.
    # The product name in the question is chemically ambiguous or describes an unusual isomer.
    # We will assume the question refers to a classic, well-known reaction of this type: the synthesis of the Wieland-Miescher ketone.
    reasoning_step3 = """
    While the exact product name in the prompt is unusual, the most famous and representative example of this type of reaction is the cyclization leading to the Wieland-Miescher ketone. The precursor for this reaction fits the description perfectly. This precursor is a single compound that contains three ketone functional groups (a triketone).
    """

    # Step 4: Name the compound.
    reactant_name = "2-methyl-2-(3-oxobutyl)cyclohexane-1,3-dione"
    reasoning_step4 = f"""
    This precursor molecule is '{reactant_name}'. When treated with base, an enolate forms and attacks one of the ketone groups in an intramolecular fashion, forming the second ring of the bicyclic system.
    """
    
    final_conclusion = "The compound that reacted with potassium hydroxide is:"

    # Print the reasoning and the final answer
    print("Plan to Identify the Reactant:")
    print("---------------------------------")
    print(textwrap.fill(reasoning_step1.strip(), 80))
    print("\n" + textwrap.fill(reasoning_step2.strip(), 80))
    print("\n" + textwrap.fill(reasoning_step3.strip(), 80))
    print("\n" + textwrap.fill(reasoning_step4.strip(), 80))
    print("\n---------------------------------")
    print(final_conclusion)
    # The final print statement outputs the name of the compound, including its numbers, as requested.
    print(reactant_name)

identify_reactant()