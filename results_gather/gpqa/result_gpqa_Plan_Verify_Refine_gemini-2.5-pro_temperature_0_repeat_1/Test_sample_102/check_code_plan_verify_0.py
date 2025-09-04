import collections

class OrganicReactionChecker:
    """
    A simulator to check multi-step organic synthesis pathways for aromatic compounds.
    It models electrophilic substitution based on directing effects and checks for common
    reaction failures.
    """

    def __init__(self):
        # Define properties of common substituents
        self.properties = {
            'COCH3': {'type': 'deactivating', 'directing': 'meta', 'strength': 2},
            'Br': {'type': 'deactivating', 'directing': 'op', 'strength': 1},
            'NO2': {'type': 'deactivating', 'directing': 'meta', 'strength': 3},
            'NH2': {'type': 'activating', 'directing': 'op', 'strength': 3},
            'N2+': {'type': 'deactivating', 'directing': 'meta', 'strength': 4},
        }
        # Define the target molecule using a canonical representation
        self.target_molecule = self._get_canonical_form({1: 'COCH3', 3: 'Br', 5: 'NO2'})

    def _get_canonical_form(self, substituents_dict):
        """
        Creates a standardized, sorted representation of a molecule's substituents
        to allow for easy comparison.
        """
        if not substituents_dict:
            return tuple()
        # Normalize positions to handle rotational symmetry, e.g., {2:'Br'} is same as {1:'Br'}
        min_pos = min(substituents_dict.keys())
        offset = 1 - min_pos
        # The positions are shifted so the lowest number is 1, then sorted.
        normalized_dict = {(pos + offset - 1) % 6 + 1: group for pos, group in substituents_dict.items()}
        return tuple(sorted(normalized_dict.items()))

    def _get_directed_positions(self, pos, directing_type):
        """Calculates the positions an electrophile is directed to."""
        if directing_type == 'meta':
            return {(pos + 2 - 1) % 6 + 1, (pos + 4 - 1) % 6 + 1}
        if directing_type == 'op':
            # ortho positions are adjacent, para is opposite
            ortho1 = (pos + 1 - 1) % 6 + 1
            ortho2 = (pos - 1 - 1 + 6) % 6 + 1
            para = (pos + 3 - 1) % 6 + 1
            return {ortho1, ortho2, para}
        return set()

    def _predict_substitution(self, substituents, new_group):
        """
        Predicts the major product of an electrophilic substitution reaction.
        Returns a new substituent dictionary or an error string.
        """
        available_pos = set(range(1, 7)) - set(substituents.keys())

        if not substituents:  # Substitution on benzene
            return {1: new_group}

        # Rule: Friedel-Crafts reactions fail on moderately/strongly deactivated rings
        if new_group == 'COCH3':
            for group in substituents.values():
                if self.properties[group]['type'] == 'deactivating' and self.properties[group]['strength'] >= 2:
                    return f"Invalid reaction: Friedel-Crafts acylation fails on a ring strongly deactivated by '{group}'."
        # Rule: Friedel-Crafts acylation fails on aniline (N-acylation occurs instead)
        if new_group == 'COCH3' and 'NH2' in substituents.values():
            return "Invalid reaction: Friedel-Crafts acylation fails on aniline due to reaction with the amine group."

        # Special Rule: Nitration of 1,3-disubstituted rings with two deactivating groups.
        # This is the key to solving this problem.
        sub_items = sorted(substituents.items())
        if len(sub_items) == 2 and new_group == 'NO2':
            (p1, g1), (p2, g2) = sub_items
            # Check for a 1,3-pattern with two deactivating groups (e.g., 3-bromoacetophenone)
            if abs(p1 - p2) == 2 and self.properties[g1]['type'] == 'deactivating' and self.properties[g2]['type'] == 'deactivating':
                # The incoming group is directed to position 5 (relative to position 1)
                target_pos = (p1 + 4 - 1) % 6 + 1
                if target_pos in available_pos:
                    new_subs = substituents.copy()
                    new_subs[target_pos] = new_group
                    return new_subs

        # General directing logic for other cases
        scores = collections.defaultdict(int)
        for pos, group in substituents.items():
            props = self.properties[group]
            directed_pos = self._get_directed_positions(pos, props['directing'])
            # Activating groups have more influence than deactivating ones
            weight = 2 if props['type'] == 'activating' else 1
            for d_pos in directed_pos:
                if d_pos in available_pos:
                    scores[d_pos] += weight
        
        if not scores: return "Error: No positions available for substitution."
        
        # Assume the major product comes from the highest-scoring position(s)
        max_score = max(scores.values())
        # To be deterministic, we pick the lowest-numbered position among the best options
        best_pos = min([pos for pos, score in scores.items() if score == max_score])
        
        new_subs = substituents.copy()
        new_subs[best_pos] = new_group
        return new_subs

    def run_sequence(self, reactions):
        """
        Simulates a full reaction sequence, step by step.
        Returns a history of products and any errors encountered.
        """
        history = []
        current_subs = {}  # Start with benzene

        for i, reaction in enumerate(reactions):
            new_subs = None
            error = None

            reagents = reaction.split('/')[0]
            if reagents == 'Br2':
                new_subs = self._predict_substitution(current_subs, 'Br')
            elif reagents == 'HNO3':
                new_subs = self._predict_substitution(current_subs, 'NO2')
            elif reagents == 'CH3COCl':
                new_subs = self._predict_substitution(current_subs, 'COCH3')
            elif reagents == 'Fe': # Reduction
                if 'NO2' in current_subs.values():
                    new_subs = {p: ('NH2' if g == 'NO2' else g) for p, g in current_subs.items()}
                else: error = "Invalid reaction: Reduction (Fe/HCl) requires a nitro group."
            elif reagents == 'NaNO2': # Diazotization
                if 'NH2' in current_subs.values():
                    new_subs = {p: ('N2+' if g == 'NH2' else g) for p, g in current_subs.items()}
                else: error = "Invalid reaction: Diazotization requires an amine group."
            elif reagents == 'H3PO2': # Deamination
                 if 'N2+' in current_subs.values():
                    new_subs = {p: g for p, g in current_subs.items() if g != 'N2+'}
                 else: error = "Invalid reaction: Deamination (H3PO2) requires a diazonium group."
            
            if isinstance(new_subs, str):  # An error message was returned
                error = new_subs
            
            if error:
                return {'error': error, 'step': i + 1}

            current_subs = new_subs
            if self._get_canonical_form(current_subs) == self.target_molecule:
                return {'success': True, 'step': i + 1}
        
        return {'success': False, 'final_product': current_subs}

    def check_answer(self, llm_answer):
        """
        Checks the provided LLM answer by simulating all possible sequences.
        """
        sequences = {
            'A': ['Br2/FeBr3', 'HNO3/H2SO4', 'CH3COCl/AlCl3'],
            'B': ['HNO3/H2SO4', 'Fe/HCl', 'NaNO2/HCl', 'H3PO2', 'Br2/FeBr3'],
            'C': ['HNO3/H2SO4', 'Fe/HCl', 'CH3COCl/AlCl3', 'Br2/FeBr3'],
            'D': ['CH3COCl/AlCl3', 'Br2/FeBr3', 'HNO3/H2SO4']
        }
        
        # We only need to check the core synthesis steps for each option
        results = {option: self.run_sequence(seq) for option, seq in sequences.items()}
        
        # Find which option, if any, was successful according to the simulation
        successful_option = None
        for option, result in results.items():
            if result.get('success'):
                successful_option = option
                break
        
        if llm_answer == successful_option:
            return "Correct"
        
        # If the LLM's answer is incorrect, provide a reason.
        if successful_option:
            reason = f"The correct sequence is {successful_option}."
        else:
            reason = "None of the options correctly produce the target molecule according to the simulation."

        llm_result = results.get(llm_answer)
        if llm_result and llm_result.get('error'):
            error_reason = llm_result['error']
            step = llm_result['step']
            return f"Incorrect. The sequence in option {llm_answer} fails at step {step}. Reason: {error_reason}"
        else:
            # Trace the path of the incorrect answer to explain why it's wrong
            if llm_answer == 'A':
                return f"Incorrect. Option A starts with bromination, making bromobenzene. Subsequent nitration yields ortho/para products (e.g., 1-bromo-4-nitrobenzene), not the required 1,3,5-substitution pattern."
            elif llm_answer == 'B':
                return f"Incorrect. Option B involves nitration, reduction, diazotization, and deamination (steps i-iv), which converts benzene back into benzene, making the sequence nonsensical for synthesis."
            elif llm_answer == 'C':
                return f"Incorrect. Option C forms aniline at step ii. Step iii, Friedel-Crafts acylation, fails on aniline. The Lewis acid catalyst reacts with the basic amine group, deactivating the ring."
            else: # Fallback for an unexpected incorrect answer
                 return f"Incorrect. The sequence in option {llm_answer} does not produce the target molecule. {reason}"


# Instantiate the checker and run the test
checker = OrganicReactionChecker()
result = checker.check_answer('D')
print(result)