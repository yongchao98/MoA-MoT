import re

def identify_reactant():
  """
  This function identifies the missing reactant in the chemical synthesis,
  prints its name, formula, and the numbers present in the formula.
  
  The overall reaction is a multi-step synthesis that converts an 
  alpha,beta-unsaturated ketone into a 5-substituted cyclohexane-1,3-dione.
  This well-known pathway involves a Michael addition of a malonic ester,
  followed by an intramolecular Claisen condensation, saponification, 
  and finally decarboxylation. The required reactant is therefore a malonic ester.
  """
  reactant_name = "Diethyl malonate"
  # The formula is often written as CH2(COOEt)2 or CH2(COOC2H5)2
  reactant_formula = "CH2(COOC2H5)2"

  print(f"The name of the missing reactant is: {reactant_name}")
  print(f"The chemical formula for this reactant is: {reactant_formula}")

  # As requested, printing each number found in the chemical formula.
  numbers = re.findall(r'\d', reactant_formula)
  
  print("\nThe numbers from the reactant's formula are:")
  for number in numbers:
      print(number)

identify_reactant()