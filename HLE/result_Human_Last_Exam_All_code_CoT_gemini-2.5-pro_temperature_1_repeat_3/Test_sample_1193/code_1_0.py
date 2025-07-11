import sys
import io

# Define a class to encapsulate the diagnostic logic
class MedicalDiagnoser:
    def __init__(self, patient_data, diagnoses):
        self.patient_data = patient_data
        self.diagnoses = diagnoses
        self.scores = {}

    def evaluate(self):
        """
        Calculates a likelihood score for each diagnosis based on clinical findings.
        """
        for key, name in self.diagnoses.items():
            score = 0
            # 1. Timeline evaluation (most important differentiator)
            if key == 'D' and self.patient_data["timeline_days"] > 7:
                # Sepsis can develop weeks post-op
                score += 50
            elif key in ['A', 'B'] and self.patient_data["timeline_days"] > 1:
                # Acute reactions are ruled out by the long timeline
                score -= 100

            # 2. Lung sounds evaluation
            if self.patient_data["lung_sounds"] == "bilateral crackles":
                if key == 'D': # Classic for ARDS from Sepsis
                    score += 40
                elif key == 'E': # Possible for cardiogenic edema
                    score += 15
                elif key in ['F', 'G', 'H']: # Not typical
                    score -= 20

            # 3. Severity of hypoxemia
            if self.patient_data["hypoxemia_o2_sat"] < 85:
                 if key == 'D': # Sepsis/ARDS causes severe hypoxemia
                    score += 30
                 elif key in ['F', 'G', 'H']: # Deconditioning does not cause this
                    score -= 30

            # 4. Penalty for vague terms
            if key in ['C', 'G']:
                score -= 50

            self.scores[key] = score

    def print_results(self):
        """
        Prints the diagnostic reasoning and the final conclusion.
        """
        # Redirect stdout to capture the print output for the "equation"
        old_stdout = sys.stdout
        sys.stdout = captured_output = io.StringIO()
        
        print("Likelihood Score Equation Components:")
        
        sorted_scores = sorted(self.scores.items(), key=lambda item: item[1], reverse=True)
        
        equation_parts = []
        for key, score in sorted_scores:
            print(f"Diagnosis '{self.diagnoses[key]}' ({key}): Score = {score}")
            equation_parts.append(f"{score}({key})")
        
        # Restore stdout
        sys.stdout = old_stdout
        
        # Get the captured output to print it
        equation_output = captured_output.getvalue()
        print(equation_output)

        final_equation = " > ".join(equation_parts)
        print("\nFinal Equation (Scores from most to least likely):")
        print(final_equation)
        
        best_choice_key = sorted_scores[0][0]
        print(f"\nConclusion: The clinical evidence most strongly supports '{self.diagnoses[best_choice_key]}'.")
        print("The 29-day timeline after major surgery, severe hypoxemia, and bilateral crackles are classic for Sepsis leading to Acute Respiratory Distress Syndrome (ARDS).")
        
        # Final answer format
        print(f"\n<<<{best_choice_key}>>>")


# Patient's clinical data
patient_vignette = {
    "timeline_days": 29,
    "hypoxemia_o2_sat": 82,
    "lung_sounds": "bilateral crackles",
    "context": "Post-Whipple procedure"
}

# Answer choices
answer_choices = {
    "A": "Acute blood transfusion reaction",
    "B": "Iodine-related reaction",
    "C": "Sensitivity reaction",
    "D": "Sepsis",
    "E": "Myocyte necrosis",
    "F": "Respiratory deconditioning",
    "G": "Lung exhaustion",
    "H": "Air pollution sensitivity"
}

# Run the diagnosis
diagnoser = MedicalDiagnoser(patient_vignette, answer_choices)
diagnoser.evaluate()
diagnoser.print_results()