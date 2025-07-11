import pandas as pd

def analyze_drug_data():
    """
    Analyzes experimental data for an antibody-drug conjugate (ADC)
    and evaluates several conclusions.
    """

    # --- Data Organization ---
    # Experiment 1: Ear Swelling (Efficacy) - Lower is better
    exp1_data = {
        'Drug': ['Anti-TNF-GRM', 'Anti-TNF-GRM', 'Anti-TNF-GRM', 'Anti-TNF-GRM',
                 'Anti-TNF', 'Anti-TNF', 'Anti-TNF', 'Anti-TNF'],
        'Dose (mg/kg)': [0.1, 1, 10, 100, 0.1, 1, 10, 100],
        'Ear Swelling (mm)': [0.04, 0.03, 0.02, 0.0, 0.4, 0.4, 0.30, 0.02]
    }
    df_exp1 = pd.DataFrame(exp1_data)

    # Experiment 2: Paw Swelling (Efficacy) - Lower/Negative is better
    exp2_data = {
        'Drug': ['Anti-TNF-GRM', 'Anti-TNF', 'GRM', 'Placebo'],
        'Day 14 Paw Swelling (mm)': [0.0, 0.5, -0.01, 0.8]
    }
    df_exp2 = pd.DataFrame(exp2_data)

    # Experiment 3: Bone Density Change (Side Effect) - More negative is worse
    exp3_data = {
        'Drug': ['Anti-TNF-GRM', 'Anti-TNF', 'GRM', 'Placebo'],
        'Dosage (mg/kg)': [10, 10, 3, 'N/A'],
        'Day 14 Bone Loss (mm^3)': [-0.3, -0.75, -0.2, -0.1]
    }
    df_exp3 = pd.DataFrame(exp3_data)

    print("--- Analysis of Experimental Data ---\n")

    # --- Evaluation of Answer Choices ---

    print("--- Evaluating Choice A: 'The ADC is less efficient...' ---")
    adc_eff_exp1 = df_exp1[(df_exp1['Drug'] == 'Anti-TNF-GRM') & (df_exp1['Dose (mg/kg)'] == 10)]['Ear Swelling (mm)'].iloc[0]
    anti_tnf_eff_exp1 = df_exp1[(df_exp1['Drug'] == 'Anti-TNF') & (df_exp1['Dose (mg/kg)'] == 10)]['Ear Swelling (mm)'].iloc[0]
    print(f"Exp 1 (10mg/kg): ADC ear swelling = {adc_eff_exp1}mm, Anti-TNF ear swelling = {anti_tnf_eff_exp1}mm.")
    adc_eff_exp2 = df_exp2[df_exp2['Drug'] == 'Anti-TNF-GRM']['Day 14 Paw Swelling (mm)'].iloc[0]
    anti_tnf_eff_exp2 = df_exp2[df_exp2['Drug'] == 'Anti-TNF']['Day 14 Paw Swelling (mm)'].iloc[0]
    print(f"Exp 2 (10mg/kg): ADC paw swelling = {adc_eff_exp2}mm, Anti-TNF paw swelling = {anti_tnf_eff_exp2}mm.")
    print("Conclusion: The ADC is significantly MORE efficient than Anti-TNF. Statement A is FALSE.\n")

    print("--- Evaluating Choices B, D, H: '...at the same risk of osteoporosis...' ---")
    adc_bone_loss = df_exp3[df_exp3['Drug'] == 'Anti-TNF-GRM']['Day 14 Bone Loss (mm^3)'].iloc[0]
    anti_tnf_bone_loss = df_exp3[df_exp3['Drug'] == 'Anti-TNF']['Day 14 Bone Loss (mm^3)'].iloc[0]
    print(f"Exp 3 (10mg/kg): ADC bone loss = {adc_bone_loss} mm^3, Anti-TNF bone loss = {anti_tnf_bone_loss} mm^3.")
    print(f"Conclusion: The bone loss from Anti-TNF ({anti_tnf_bone_loss}) is 2.5 times greater than from ADC ({adc_bone_loss}). The risks are not the same. Statements B, D, and H are FALSE.\n")

    print("--- Evaluating Choice G: '...ADC but not GMR can fight inflamaiton.' ---")
    grm_eff_exp2 = df_exp2[df_exp2['Drug'] == 'GRM']['Day 14 Paw Swelling (mm)'].iloc[0]
    placebo_eff_exp2 = df_exp2[df_exp2['Drug'] == 'Placebo']['Day 14 Paw Swelling (mm)'].iloc[0]
    print(f"Exp 2: GRM resulted in paw swelling of {grm_eff_exp2}mm, compared to Placebo at {placebo_eff_exp2}mm.")
    print("Conclusion: GRM is effective against inflammation. Statement G is FALSE.\n")

    print("--- Evaluating Choices F, I: 'GRM will induce fewer side effects than... ADC... at the same dosage.' ---")
    grm_bone_loss = df_exp3[df_exp3['Drug'] == 'GRM']['Day 14 Bone Loss (mm^3)'].iloc[0]
    grm_dose = df_exp3[df_exp3['Drug'] == 'GRM']['Dosage (mg/kg)'].iloc[0]
    adc_bone_loss = df_exp3[df_exp3['Drug'] == 'Anti-TNF-GRM']['Day 14 Bone Loss (mm^3)'].iloc[0]
    adc_dose = df_exp3[df_exp3['Drug'] == 'Anti-TNF-GRM']['Dosage (mg/kg)'].iloc[0]
    grm_side_effect_per_mg = grm_bone_loss / grm_dose
    adc_side_effect_per_mg = adc_bone_loss / adc_dose
    print(f"Side effect per mg/kg for GRM: {grm_bone_loss}/{grm_dose} = {grm_side_effect_per_mg:.3f}")
    print(f"Side effect per mg/kg for ADC: {adc_bone_loss}/{adc_dose} = {adc_side_effect_per_mg:.3f}")
    print("Conclusion: GRM has more potent side effects per mg than the ADC. Extrapolating to the same 10mg/kg dose suggests GRM would have WORSE side effects. Statements F and I are FALSE.\n")

    print("--- Evaluating Choice E: 'The dosage... was chosen correctly to compare...' ---")
    adc_equipotent = df_exp1[(df_exp1['Drug'] == 'Anti-TNF-GRM') & (df_exp1['Dose (mg/kg)'] == 10)]
    anti_tnf_equipotent = df_exp1[(df_exp1['Drug'] == 'Anti-TNF') & (df_exp1['Dose (mg/kg)'] == 100)]
    print(f"Exp 1 shows ADC at {adc_equipotent['Dose (mg/kg)'].iloc[0]} mg/kg has an effect of {adc_equipotent['Ear Swelling (mm)'].iloc[0]}mm.")
    print(f"Exp 1 shows Anti-TNF requires {anti_tnf_equipotent['Dose (mg/kg)'].iloc[0]} mg/kg for the same effect ({anti_tnf_equipotent['Ear Swelling (mm)'].iloc[0]}mm).")
    print("Conclusion: The doses used for side-effect comparison (10mg/kg) are not equipotent. A rigorous comparison would use equipotent doses. Therefore, the choice of dosage can be considered scientifically suboptimal or 'incorrect'. Statement E is likely FALSE.\n")
    
    print("--- Final Conclusion ---")
    print("Based on the analysis, statements A, B, D, E, F, G, H, and I are all determined to be false or incorrect.")
    print("Therefore, the correct option is C.")

if __name__ == '__main__':
    analyze_drug_data()
<<<C>>>