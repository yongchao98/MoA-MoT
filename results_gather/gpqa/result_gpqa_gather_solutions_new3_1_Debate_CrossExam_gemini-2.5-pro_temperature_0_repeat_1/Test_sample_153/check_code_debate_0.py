import sys
from io import StringIO

class ChemistryChecker:
    """
    A class to check the correctness of a proposed chemical structure based on spectral data.
    """
    def __init__(self, proposed_answer_key):
        self.compounds = {
            "A": {
                "name": "4-chlorobenzoic acid",
                "formula": "C7H5ClO2",
                "functional_groups": ["carboxylic acid", "aromatic ring", "halide"],
                "substitution": "para",  # 1,4-disubstituted
                "aromatic_protons": 4
            },
            "B": {
                "name": "Phenyl chloroformate",
                "formula": "C7H5ClO2",
                "functional_groups": ["chloroformate", "aromatic ring", "halide"],
                "substitution": "mono",
                "aromatic_protons": 5
            },
            "C": {
                "name": "3-Chloro-2-hydroxybenzaldehyde",
                "formula": "C7H5ClO2",
                "functional_groups": ["aldehyde", "phenol", "aromatic ring", "halide"],
                "substitution": "1,2,3-trisubstituted",
                "aromatic_protons": 3
            },
            "D": {
                "name": "2-chlorobenzoic acid",
                "formula": "C7H5ClO2",
                "functional_groups": ["carboxylic acid", "aromatic ring", "halide"],
                "substitution": "ortho",  # 1,2-disubstituted
                "aromatic_protons": 4
            }
        }
        # The final answer from the LLM to be checked
        self.proposed_answer_key = proposed_answer_key
        self.errors = []

    def calculate_mw_35Cl(self, formula):
        """Calculates molecular weight using the 35Cl isotope."""
        atomic_masses_35Cl = {'C': 12.000, 'H': 1.008, 'O': 15.999, 'Cl': 34.969}
        # Formula is C7H5ClO2 for all options
        mw = (7 * atomic_masses_35Cl['C'] +
              5 * atomic_masses_35Cl['H'] +
              1 * atomic_masses_35Cl['Cl'] +
              2 * atomic_masses_35Cl['O'])
        return mw

    def check_ms(self):
        """Checks consistency with Mass Spectrometry data."""
        target_mw = 156
        compound_data = self.compounds[self.proposed_answer_key]
        
        # Check 1: Molecular Weight
        calculated_mw = self.calculate_mw_35Cl(compound_data["formula"])
        if not (target_mw - 1 < calculated_mw < target_mw + 1):
            self.errors.append(f"MS Check Failed: Calculated molecular weight for {compound_data['name']} ({calculated_mw:.2f}) does not match the M+ peak at m/z = {target_mw}.")

        # Check 2: Presence of Chlorine for M+2 peak
        if "halide" not in compound_data["functional_groups"] or "Cl" not in compound_data["formula"]:
             self.errors.append(f"MS Check Failed: The M+2 peak at m/z=158 (32% intensity) indicates a Chlorine atom, but the proposed structure {compound_data['name']} does not contain Chlorine.")

    def check_ir(self):
        """Checks consistency with Infrared Spectroscopy data."""
        # The combination of a very broad peak (3500-2700 cm-1) and a C=O peak (1720 cm-1)
        # is characteristic of a carboxylic acid.
        compound_data = self.compounds[self.proposed_answer_key]
        
        if "carboxylic acid" not in compound_data["functional_groups"]:
            self.errors.append(f"IR Check Failed: The spectrum is characteristic of a carboxylic acid. The proposed structure, {compound_data['name']}, is not a carboxylic acid.")

    def check_nmr(self):
        """Checks consistency with 1H NMR Spectroscopy data."""
        compound_data = self.compounds[self.proposed_answer_key]
        
        # Check 1: Carboxylic acid proton signal at 11.0 ppm
        if "carboxylic acid" not in compound_data["functional_groups"]:
            self.errors.append(f"NMR Check Failed: The signal at 11.0 ppm is characteristic of a carboxylic acid proton, but the proposed structure, {compound_data['name']}, is not a carboxylic acid.")
            
        # Check 2: Aromatic protons count
        total_aromatic_protons_data = 4 # from 2H + 2H
        if compound_data.get("aromatic_protons") != total_aromatic_protons_data:
            self.errors.append(f"NMR Check Failed: The data shows {total_aromatic_protons_data} aromatic protons, but the proposed structure, {compound_data['name']}, has {compound_data['aromatic_protons']}.")

        # Check 3: Para-substitution pattern
        # The pattern of two doublets (2H each) is a classic signature for a para-substituted ring.
        if compound_data.get("substitution") != "para":
            self.errors.append(f"NMR Check Failed: The aromatic region shows a pattern for a para-substituted ring. The proposed structure, {compound_data['name']}, has '{compound_data['substitution']}' substitution, which would produce a different splitting pattern.")

    def run_checks(self):
        """Runs all checks and prints the final result."""
        self.check_ms()
        self.check_ir()
        self.check_nmr()
        
        if not self.errors:
            print("Correct")
        else:
            print("Incorrect")
            for error in self.errors:
                print(f"- {error}")

# The final proposed answer is 'A'.
# We instantiate the checker with 'A' and run it.
try:
    checker = ChemistryChecker(proposed_answer_key="A")
    checker.run_checks()
except KeyError:
    print("Incorrect: The proposed answer key is not a valid option (A, B, C, or D).")
except Exception as e:
    print(f"An error occurred during the check: {e}")
