def solve_nmr_puzzle():
  """
  This function determines the splitting pattern and integration of the most
  deshielded proton in Compound 1.

  The reasoning is as follows:
  1. Reaction: The reaction is an electrophilic aromatic sulfonation of the
     di-propyl diazaoxatriangulenium (Pr-DAOTA) cation with concentrated sulfuric acid.
     This introduces sulfonic acid (-SO3H) groups onto the aromatic rings.
  2. Product Structure: Sulfonation occurs at the positions para to the activating
     nitrogen atoms, leading to a symmetrical disulfonated product.
  3. Most Deshielded Proton: The proton on the central aromatic ring, positioned
     between the two electron-withdrawing nitrogen atoms in the cationic core, is the
     most electron-deficient and therefore the most deshielded proton.
  4. Splitting Pattern: This proton has no protons on adjacent carbons, so it does
     not experience spin-spin coupling. Thus, its signal is a singlet.
  5. Integration: There is only one proton of this type in the molecule. Therefore,
     the integration of its peak corresponds to one proton (1H).
  """
  splitting_pattern = "singlet"
  integration_value = 1
  integration_unit = "H"

  print(f"The highest deshielded proton peak in Compound 1 is a {splitting_pattern}.")
  print(f"The integration of this peak is {integration_value}{integration_unit}.")

solve_nmr_puzzle()