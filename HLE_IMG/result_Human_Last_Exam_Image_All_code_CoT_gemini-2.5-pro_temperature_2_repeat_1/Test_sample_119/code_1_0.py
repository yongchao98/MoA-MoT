def solve_spectroscopy_problem():
    """
    Analyzes the provided spectroscopic data to determine the IUPAC name of the compound.
    """
    analysis_steps = """
    Step-by-step analysis to determine the IUPAC name:

    1. Molecular Formula and Degree of Unsaturation (DoU)
    - Mass Spectrum (MS): The molecular ion peak (M+) at m/z 135 indicates a molecular weight of 135 g/mol.
    - Nitrogen Rule: An odd molecular weight suggests the presence of an odd number of nitrogen atoms. Assuming one nitrogen (N).
    - Molecular Formula Calculation:
        - Mass remaining after N: 135 - 14 = 121.
        - Number of carbons (C): 121 / 12 is ~10. Let's try 9 carbons. 9 * 12 = 108.
        - Number of hydrogens (H): 121 - 108 = 13.
        - The molecular formula is C9H13N.
    - Degree of Unsaturation (DoU):
        - DoU = C - H/2 + N/2 + 1
        - DoU = 9 - (13/2) + (1/2) + 1 = 9 - 6.5 + 0.5 + 1 = 4.
        - A DoU of 4 is characteristic of a monosubstituted benzene ring.

    2. Functional Group Analysis (IR and NMR)
    - IR Spectrum: Shows key stretches for N-H (primary amine, ~3300 cm⁻¹), aromatic C-H (>3000 cm⁻¹), and aliphatic C-H (<3000 cm⁻¹).
    - ¹³C NMR / DEPT-135: 7 total carbon signals.
        - DEPT-135 shows one negative signal (one CH2 group) and five positive signals (CH and/or CH3 groups).
        - δ 145.1: Quaternary C (aromatic ipso, not in DEPT).
        - δ 128.5, 127.3, 126.3: Aromatic CH groups (positive).
        - δ 49.6: Aliphatic CH group (positive).
        - δ 43.5: Aliphatic CH2 group (negative).
        - δ 19.2: Aliphatic CH3 group (positive).
        - This confirms a C9 skeleton: 1 C(quat), 5 CH(arom), 1 CH(aliph), 1 CH2(aliph), 1 CH3(aliph).

    3. Structure Assembly from 1H, 13C, and HSQC data
    - The data points to these fragments: a C6H5- group, an -NH2 group, and a C3 aliphatic chain.
    - The aliphatic chain consists of one CH, one CH2, and one CH3.
    - 1H NMR signals confirm the fragments:
        - δ ~7.2 ppm (5H): A monosubstituted phenyl group (C6H5-).
        - δ ~1.2 ppm (3H, doublet): A -CH3 group next to a CH group.
        - δ ~2.7-2.9 ppm (3H total): A -CH- group (1H) and a -CH2- group (2H) that are adjacent.
    - HSQC confirms the connectivity:
        - H(~1.2 ppm) connects to C(19.2 ppm) -> This is the CH3 group.
        - H(~2.7 ppm) connects to C(43.5 ppm) -> This is the CH2 group.
        - H(~2.8 ppm) connects to C(49.6 ppm) -> This is the CH group.
    - Putting the fragments together, the structure must be C6H5-CH2-CH(NH2)-CH3.
    - This structure (1-phenylpropan-2-amine) perfectly matches all spectral data:
        - The benzylic CH2 protons at ~2.7 ppm are attached to the CH2 carbon at 43.5 ppm.
        - The CH proton (alpha to the amine) at ~2.8 ppm is attached to the CH carbon at 49.6 ppm.
        - The CH3 protons at ~1.2 ppm are attached to the CH3 carbon at 19.2 ppm and appear as a doublet due to splitting by the neighboring single CH proton.

    4. IUPAC Name
    - The parent chain is propane.
    - Numbering gives the principal functional group (amine) the lowest number (position 2). The phenyl substituent is at position 1.
    - The final IUPAC name is 1-phenylpropan-2-amine.
    """
    print(analysis_steps)

solve_spectroscopy_problem()
<<<1-phenylpropan-2-amine>>>