def compare_toxicity():
    """
    Compares the acute toxicity of Trimethyltin chloride (TMT-Cl) and
    Tributyltin chloride (TBT-Cl) using their LD50 values.
    """
    # LD50 values are a standard measure of acute toxicity.
    # The data represents the oral LD50 in rats (mg of substance per kg of body weight).
    # A lower LD50 value signifies higher toxicity.
    tmt_cl_ld50 = 12.6  # mg/kg
    tbt_cl_ld50 = 129.0  # mg/kg

    print("To understand the relative danger of chemical compounds, we often use the LD50 value.")
    print("LD50 (Lethal Dose, 50%) is the dose that is lethal to 50% of a test population.")
    print("A lower LD50 value indicates a higher toxicity.\n")

    print("--- Toxicity Comparison ---")
    print(f"Trimethyltin chloride (TMT-Cl) LD50: {tmt_cl_ld50} mg/kg")
    print(f"Tributyltin chloride (TBT-Cl) LD50: {tbt_cl_ld50} mg/kg")
    print("---------------------------\n")

    # Conclusion based on the data
    print("Conclusion:")
    print(f"The LD50 of Trimethyltin chloride (TMT-Cl) is {tmt_cl_ld50} mg/kg, which is significantly lower than the LD50 of Tributyltin chloride (TBT-Cl) at {tbt_cl_ld50} mg/kg.")
    print("This means a much smaller amount of TMT-Cl is required to cause a lethal effect, making it far more acutely dangerous.")
    print("\nTherefore, the most important factor demonstrating TMT-Cl's higher danger is its significantly lower LD50 value.")

if __name__ == '__main__':
    compare_toxicity()