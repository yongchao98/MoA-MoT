import re

class ChemicalProblemChecker:
    """
    This class checks the correctness of the provided answer by verifying
    each constraint given in the chemistry problem statement.
    """

    # A small database of properties for the compounds involved in the problem.
    # This encapsulates the chemical knowledge required to solve the problem.
    COMPOUNDS_DB = {
        'cyclohexane': {
            'formula': 'C6H12',
            'is_saturated': True,
            'is_aromatic': False,
            'has_conjugated_bonds': False,
        },
        'benzene': {
            'formula': 'C6H6',
            'is_saturated': False,
            'is_aromatic': True,
            'has_conjugated_bonds': True,  # Aromaticity is a special case of conjugation
        },
        'cyclohexene': {
            'formula': 'C6H10',
            'is_saturated': False,
            'is_aromatic': False,
            'has_conjugated_bonds': False,
        },
        '1,3-cyclohexadiene': {
            'formula': 'C6H8',
            'is_saturated': False,
            'is_aromatic': False,
            'has_conjugated_bonds': True,
        },
        '1,4-cyclohexadiene': {
            'formula': 'C6H8',
            'is_saturated': False,
            'is_aromatic': False,
            'has_conjugated_bonds': False,
        }
    }
    ATOMIC_MASSES = {'C': 12.011, 'H': 1.008}

    def _get_atoms(self, formula: str) -> dict:
        """Parses a formula like 'C6H12' into a dict {'C': 6, 'H': 12}."""
        match = re.match(r'C(\d+)H(\d+)', formula)
        if not match:
            raise ValueError(f"Invalid formula format: {formula}")
        return {'C': int(match.group(1)), 'H': int(match.group(2))}

    def _calculate_h_mass_fraction(self, formula: str) -> float:
        """Calculates the mass fraction of hydrogen in a given hydrocarbon formula."""
        atoms = self._get_atoms(formula)
        mass_h = atoms['H'] * self.ATOMIC_MASSES['H']
        mass_c = atoms['C'] * self.ATOMIC_MASSES['C']
        if (mass_c + mass_h) == 0: return 0.0
        return mass_h / (mass_c + mass_h)

    def _decolorizes_bromine_water(self, compound_name: str) -> bool:
        """Checks if a compound decolorizes bromine water. True for unsaturated, non-aromatic."""
        props = self.COMPOUNDS_DB[compound_name]
        return not props['is_saturated'] and not props['is_aromatic']

    def _hydrogenation_product(self, compound_name: str) -> str:
        """Determines the hydrogenation product. For this problem, all C6 rings go to cyclohexane."""
        if 'C6' in self.COMPOUNDS_DB[compound_name]['formula']:
            return 'cyclohexane'
        return None

    def check_answer(self):
        """
        This method runs a full check on the proposed solution based on the problem's constraints.
        The proposed solution is derived from the provided LLM answer.
        """
        # Based on the provided answer's reasoning:
        proposed_Z_name = 'cyclohexane'
        proposed_Y_components = ['cyclohexane', 'benzene']
        proposed_X_components = ['cyclohexene', '1,4-cyclohexadiene']
        # The final answer is 18 (Option C)
        final_answer_value = 18

        errors = []

        # --- Step 1: Verify Substance Z ---
        z_props = self.COMPOUNDS_DB[proposed_Z_name]
        # Constraint: Mass fraction of hydrogen is 14.28%
        h_frac = self._calculate_h_mass_fraction(z_props['formula'])
        if not (0.142 <= h_frac <= 0.144):  # Check if it's ~14.3%
            errors.append(f"Constraint Fail (Z): {proposed_Z_name} has H mass fraction {h_frac:.2%}, which is not ~14.28%.")
        # Constraint: Z does not react further with hydrogen (is saturated)
        if not z_props['is_saturated']:
            errors.append(f"Constraint Fail (Z): {proposed_Z_name} is not saturated, but the problem states it is.")

        # --- Step 2: Verify Mixture Y ---
        # Constraint: Y does not decolorize bromine water
        for comp in proposed_Y_components:
            if self._decolorizes_bromine_water(comp):
                errors.append(f"Constraint Fail (Y): Component '{comp}' of mixture Y should not decolorize bromine water, but it does.")
        # Constraint: Hydrogenation of Y gives only Z
        for comp in proposed_Y_components:
            if self._hydrogenation_product(comp) != proposed_Z_name:
                errors.append(f"Constraint Fail (Y): Hydrogenation of '{comp}' from Y does not yield Z ({proposed_Z_name}).")
        # Constraint: Z is a constituent of Y
        if proposed_Z_name not in proposed_Y_components:
            errors.append(f"Constraint Fail (Y): Z ({proposed_Z_name}) is not a constituent of mixture Y.")

        # --- Step 3: Verify Mixture X ---
        # Constraint: X decolorizes bromine water
        for comp in proposed_X_components:
            if not self._decolorizes_bromine_water(comp):
                errors.append(f"Constraint Fail (X): Component '{comp}' of mixture X should decolorize bromine water, but it does not.")
        # Constraint: No conjugated multiple bonds in X
        for comp in proposed_X_components:
            if self.COMPOUNDS_DB[comp]['has_conjugated_bonds']:
                errors.append(f"Constraint Fail (X): Component '{comp}' of mixture X has conjugated bonds, which is forbidden.")
        # Constraint: Hydrogenation of X gives only Z
        for comp in proposed_X_components:
            if self._hydrogenation_product(comp) != proposed_Z_name:
                errors.append(f"Constraint Fail (X): Hydrogenation of '{comp}' from X does not yield Z ({proposed_Z_name}).")

        # --- Step 4: Verify the Disproportionation Reaction (Atom Conservation) ---
        x1_atoms = self._get_atoms(self.COMPOUNDS_DB[proposed_X_components[0]]['formula'])
        x2_atoms = self._get_atoms(self.COMPOUNDS_DB[proposed_X_components[1]]['formula'])
        y1_atoms = self._get_atoms(self.COMPOUNDS_DB[proposed_Y_components[0]]['formula'])
        y2_atoms = self._get_atoms(self.COMPOUNDS_DB[proposed_Y_components[1]]['formula'])

        total_x_h = x1_atoms['H'] + x2_atoms['H']
        total_y_h = y1_atoms['H'] + y2_atoms['H']
        total_x_c = x1_atoms['C'] + x2_atoms['C']
        total_y_c = y1_atoms['C'] + y2_atoms['C']

        if total_x_h != total_y_h or total_x_c != total_y_c:
            errors.append(f"Constraint Fail (Reaction): Atom conservation fails in X -> Y reaction. Reactants: C{total_x_c}H{total_x_h}, Products: C{total_y_c}H{total_y_h}.")

        # --- Step 5: Verify the Final Answer ---
        # The question asks for the total number of hydrogen atoms in the two liquids of mixture X.
        calculated_h_count = total_x_h
        if calculated_h_count != final_answer_value:
            errors.append(f"Final Answer Mismatch: The total H atoms in the proposed components of X is {calculated_h_count}, but the final answer is {final_answer_value}.")

        # --- Final Verdict ---
        if not errors:
            return "Correct"
        else:
            return "Incorrect. Reasons:\n" + "\n".join(errors)

# Create an instance of the checker and run the check.
checker = ChemicalProblemChecker()
result = checker.check_answer()
print(result)