import math

def analyze_osteoporosis_risk():
    """
    Analyzes and prints the data to evaluate the statements in option F.
    """

    # Experiment 3 Data: Bone density change in cubic millimeters at 14 days
    exp3_data = {
        "Anti-TNF-GRM": {"dose_mg_per_kg": 10, "bone_loss_mm3": -0.3},
        "Anti-TNF": {"dose_mg_per_kg": 10, "bone_loss_mm3": -0.75},
        "GRM": {"dose_mg_per_kg": 3, "bone_loss_mm3": -0.2},
        "Placebo": {"dose_mg_per_kg": None, "bone_loss_mm3": -0.1}
    }

    # --- Part 1: "The mice treated with anti-TNF are at risk of osteoporosis." ---
    print("--- Analysis for Statement F ---")
    print("\nPart 1: Is anti-TNF associated with osteoporosis risk?")
    anti_tnf_loss = exp3_data["Anti-TNF"]["bone_loss_mm3"]
    placebo_loss = exp3_data["Placebo"]["bone_loss_mm3"]
    print(f"Anti-TNF bone loss: {anti_tnf_loss} mm3")
    print(f"Placebo bone loss: {placebo_loss} mm3")
    is_at_risk = anti_tnf_loss < placebo_loss
    print(f"Conclusion: Since {anti_tnf_loss} < {placebo_loss}, the statement that anti-TNF treated mice are at risk is TRUE.")

    # --- Part 2: "The side effects of the tested ADC are lower than those of the anti-TFN." ---
    print("\nPart 2: Are ADC side effects lower than anti-TNF side effects?")
    adc_loss = exp3_data["Anti-TNF-GRM"]["bone_loss_mm3"]
    adc_dose = exp3_data["Anti-TNF-GRM"]["dose_mg_per_kg"]
    anti_tnf_dose = exp3_data["Anti-TNF"]["dose_mg_per_kg"]
    
    # Check if doses are comparable
    if adc_dose == anti_tnf_dose:
        print(f"Comparing ADC and anti-TNF at the same dose of {adc_dose} mg/kg:")
        print(f"ADC (Anti-TNF-GRM) bone loss: {adc_loss} mm3")
        print(f"Anti-TNF bone loss: {anti_tnf_loss} mm3")
        side_effects_lower = adc_loss > anti_tnf_loss  # Less bone loss means a greater number
        print(f"Conclusion: Since {adc_loss} > {anti_tnf_loss}, the statement that ADC has lower side effects is TRUE.")
    else:
        print("Doses are not comparable for a direct conclusion.")


    # --- Part 3: "GRM will induce fewer side effects than the tested ADC even when the dosage will be the same." ---
    print("\nPart 3: Will GRM have fewer side effects than ADC at the same dose?")
    grm_loss = exp3_data["GRM"]["bone_loss_mm3"]
    grm_dose = exp3_data["GRM"]["dose_mg_per_kg"]
    print(f"The experiment measured GRM bone loss as {grm_loss} mm3 at a dose of {grm_dose} mg/kg.")
    print(f"It measured ADC bone loss as {adc_loss} mm3 at a dose of {adc_dose} mg/kg.")
    print("Conclusion: The data is insufficient to make a direct comparison at the same dose.")
    
    extrapolated_grm_loss = (grm_loss / grm_dose) * adc_dose
    print(f"A simple linear extrapolation for GRM at {adc_dose} mg/kg would be ({grm_loss} / {grm_dose}) * {adc_dose} = {extrapolated_grm_loss:.2f} mm3.")
    print(f"Since the extrapolated value {extrapolated_grm_loss:.2f} is more negative than the ADC's value of {adc_loss}, this statement is NOT supported by the data and is likely FALSE.")

analyze_osteoporosis_risk()
<<<F>>>