import collections

def get_reactant_info():
    """
    Identifies the reactant and provides its details based on chemical analysis.
    The reaction is a Robinson-type annulation involving a Michael addition
    of a malonic ester to an alpha,beta-unsaturated ketone, followed by
    an intramolecular Dieckmann condensation, hydrolysis, and decarboxylation
    to yield a 5-substituted cyclohexane-1,3-dione.
    """
    reactant_name = "Diethyl malonate"
    
    # Chemical formula C7H12O4
    formula = collections.OrderedDict([
        ('C', 7),
        ('H', 12),
        ('O', 4)
    ])
    
    print(f"The required reactant is: {reactant_name}")
    
    # Building and printing the chemical formula string
    formula_str = "".join([f"{atom}{count}" for atom, count in formula.items()])
    print(f"Chemical Formula: {formula_str}")

    # "output each number in the final equation!" -> interpreting this as outputting the counts in the formula.
    print("The numbers in the chemical formula are:")
    for atom, count in formula.items():
        print(f"{atom}: {count}")

get_reactant_info()