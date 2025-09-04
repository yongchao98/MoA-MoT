import sys
import io

def check_molecular_symmetry():
    """
    Checks the correctness of the LLM's answer by computationally determining
    the point group of the selected molecule.
    """
    # Suppress verbose output from libraries by redirecting stdout
    old_stdout = sys.stdout
    sys.stdout = captured_output = io.StringIO()

    try:
        import pubchempy as pcp
        from pyscf import gto
    except ImportError:
        sys.stdout = old_stdout
        return "Execution failed: Required libraries are not installed. Please run 'pip install pubchempy pyscf' to execute this check."

    # Data from the question and the provided answer
    # The LLM selected option C.
    molecule_name = "mellitic trianhydride"
    # PubChem Compound ID (CID) for mellitic trianhydride
    molecule_cid = 78261
    # The symmetry group claimed by the LLM for this molecule
    claimed_symmetry = "C3h"
    
    actual_symmetry = None
    error_message = None

    try:
        # Fetch 3D structure from PubChem
        compounds = pcp.get_compounds(molecule_cid, 'cid', record_type='3d')
        if not compounds:
            raise ValueError(f"Could not find a 3D structure for {molecule_name} (CID: {molecule_cid}) on PubChem.")
        
        compound = compounds[0]
        
        # Prepare atoms list in PySCF's format: [('AtomSymbol', (x, y, z)), ...]
        atoms_pyscf_format = []
        for atom in compound.atoms:
            atoms_pyscf_format.append((atom.element, (atom.x, atom.y, atom.z)))
        
        # Use PySCF to build a molecule object and determine its point group
        mol = gto.Mole()
        mol.atom = atoms_pyscf_format
        # The build() method automatically detects and stores the point group in mol.groupname
        mol.build(verbose=0)
        
        # Standardize the output format (e.g., 'd3h' -> 'D3h')
        actual_symmetry = mol.groupname.capitalize()

    except Exception as e:
        error_message = f"An error occurred during the verification process: {e}"

    # Restore stdout
    sys.stdout = old_stdout

    # Final evaluation
    if error_message:
        return error_message
        
    if actual_symmetry == claimed_symmetry.capitalize():
        # This case is highly unlikely based on chemical theory but is included for completeness.
        # If this were true, we would need to check the other molecules for uniqueness.
        return "Correct"
    else:
        reason = (
            f"Incorrect. The answer claims that molecule C (mellitic trianhydride) has C3h symmetry. "
            f"However, a computational check using its 3D structure from PubChem (CID: {molecule_cid}) and the PySCF library "
            f"determines its point group to be {actual_symmetry}. "
            f"A molecule with {actual_symmetry} symmetry possesses a C3 axis and a horizontal mirror plane (σh), but it also has three C2 axes "
            f"perpendicular to the principal C3 axis. The C3h group explicitly lacks these perpendicular C2 axes. "
            f"The LLM's reasoning incorrectly states that these perpendicular C2 axes are absent, leading to the wrong point group assignment."
        )
        return reason

# Run the check and print the result
result = check_molecular_symmetry()
print(result)