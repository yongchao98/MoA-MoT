class ReactionChecker:
    """
    A class to verify the starting material in a Ring-Opening Cross-Metathesis (ROCM) reaction.
    It uses a rule-based approach based on chemical principles to predict reaction outcomes.
    """
    def __init__(self):
        # Define the key structural features of the target product.
        self.target_product_features = {
            "ring": "cyclopentane",
            "substituents": sorted(["vinyl", "prop-1-en-1-yl"]),
            "substitution_pattern": "1,2-disubstituted"
        }
        # Define the provided options.
        self.options = {
            "A": "2-methylbicyclo[3.1.0]hex-2-ene",
            "B": "2-methyl-3-methylenebicyclo[2.1.0]pentane",
            "C": "bicyclo[3.2.0]hept-6-ene",
            "D": "1,2-dimethylenecyclopentane"
        }

    def _predict_outcome(self, starting_material_name):
        """
        Predicts the outcome of the reaction for a given starting material.
        Returns a dictionary of predicted product features or an error message.
        """
        if starting_material_name == "2-methylbicyclo[3.1.0]hex-2-ene":
            # Ring-opening metathesis (ROM) of the double bond in the 5-membered ring
            # would destroy the cyclopentane core, which is required in the product.
            return {"error": "The reaction would open the cyclopentane ring, which is inconsistent with the product."}

        elif starting_material_name == "2-methyl-3-methylenebicyclo[2.1.0]pentane":
            # This starting material does not contain a cyclopentane ring and is unlikely to form one
            # with the correct substituents via a simple ROCM pathway.
            return {"error": "This starting material does not have a cyclopentane ring and is unlikely to form one under these conditions."}

        elif starting_material_name == "bicyclo[3.2.0]hept-6-ene":
            # This is the correct substrate for ROCM.
            # Structure: cyclopentane fused to a strained cyclobutene ring.
            # 1. Ring-opening of the cyclobutene preserves the cyclopentane ring.
            # 2. The mechanism correctly installs a vinyl group (from the methyleneruthenium catalyst)
            #    and a prop-1-en-1-yl group (from cross-metathesis with 1-propene)
            #    on adjacent carbons of the cyclopentane ring.
            return {
                "ring": "cyclopentane",
                "substituents": sorted(["vinyl", "prop-1-en-1-yl"]),
                "substitution_pattern": "1,2-disubstituted"
            }

        elif starting_material_name == "1,2-dimethylenecyclopentane":
            # This is a diene that undergoes simple cross-metathesis, not ROCM.
            # Reacting with 1-propene gives a prop-1-en-1-yl group, but the other group remains
            # a methylene (=CH2), not a vinyl (-CH=CH2) group.
            return {
                "ring": "cyclopentane",
                "substituents": sorted(["methylene", "prop-1-en-1-yl"]),
                "substitution_pattern": "1,2-disubstituted"
            }
        
        return {"error": "Unknown starting material."}

    def check_answer(self, proposed_answer_letter):
        """
        Checks if the proposed answer letter corresponds to the correct starting material.
        """
        if proposed_answer_letter not in self.options:
            return f"Invalid option '{proposed_answer_letter}'. The valid options are {list(self.options.keys())}."

        proposed_material_name = self.options[proposed_answer_letter]
        predicted_features = self._predict_outcome(proposed_material_name)

        # Check if the predicted features exactly match the target product features
        if predicted_features == self.target_product_features:
            return "Correct"
        else:
            # The proposed answer is incorrect. Formulate a clear reason.
            if "error" in predicted_features:
                reason = predicted_features["error"]
            else:
                predicted_subs = ", ".join(predicted_features.get("substituents", []))
                target_subs = ", ".join(self.target_product_features.get("substituents", []))
                reason = f"The predicted substituents ({predicted_subs}) do not match the target substituents ({target_subs})."
            
            return (f"The answer '{proposed_answer_letter}' is incorrect. "
                    f"Starting with '{proposed_material_name}' would not yield the target product. "
                    f"Reason: {reason}")

# --- Main execution block for the check ---
try:
    # The answer from the other LLM to be checked
    llm_answer = "C"
    
    # Instantiate the checker and run the check
    checker = ReactionChecker()
    result = checker.check_answer(llm_answer)
    
    # Print the final result
    print(result)

except Exception as e:
    print(f"An error occurred during the check: {e}")