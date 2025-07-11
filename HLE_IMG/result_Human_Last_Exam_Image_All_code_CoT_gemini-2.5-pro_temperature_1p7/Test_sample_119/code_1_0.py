def solve_structure():
    """
    This function prints the determined IUPAC name based on the provided spectral data.
    The analysis steps are as follows:
    1. MS: M+ at m/z 135 and the nitrogen rule suggest the formula C9H13N.
    2. DoU: Degree of unsaturation is 4, suggesting a benzene ring.
    3. IR: Confirms aromatic C-H, aliphatic C-H, and N-H stretches.
    4. 13C & DEPT: Reveals 4 aromatic carbons (1 quaternary C, 3 CH) and 3 aliphatic carbons (1 CH, 1 CH2, 1 CH3).
    5. 1H & HSQC: Correlates protons to carbons, identifying a C6H5- group, a -CH(CH3)- fragment, and a -CH2- group.
    6. Proton Count: 11 protons are attached to carbons, the remaining 2 protons are from an -NH2 group.
    7. Assembly: Combining the fragments (C6H5-, -CH(CH3)-, and -CH2NH2) leads to the structure 2-phenylpropan-1-amine.
    """
    iupac_name = "2-phenylpropan-1-amine"
    print(f"The IUPAC name of the compound is: {iupac_name}")

solve_structure()