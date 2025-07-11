def calculate_map():
    """
    Calculates and prints the Mean Arterial Pressure (MAP) based on given vital signs.
    """
    # Vital signs provided
    systolic_bp = 90  # mmHg
    diastolic_bp = 60 # mmHg
    heart_rate = 100 # bpm
    respiratory_rate = 40 # breaths/min

    # MAP is calculated as: (Systolic BP + 2 * Diastolic BP) / 3
    map_pressure = (systolic_bp + 2 * diastolic_bp) / 3

    print("The patient's blood pressure is {}/{} mmHg.".format(systolic_bp, diastolic_bp))
    print("The formula for Mean Arterial Pressure (MAP) is: (Systolic BP + 2 * Diastolic BP) / 3")
    print("Calculating the patient's MAP:")
    print("MAP = ({} + 2 * {}) / 3".format(systolic_bp, diastolic_bp))
    print("MAP = ({} + {}) / 3".format(systolic_bp, 2 * diastolic_bp))
    print("MAP = {} / 3".format(systolic_bp + 2 * diastolic_bp))
    print("The patient's MAP is {:.2f} mmHg.".format(map_pressure))

    print("\nOther relevant vital signs:")
    print("Heart Rate: {} bpm".format(heart_rate))
    print("Respiratory Rate: {} breaths/min".format(respiratory_rate))
    
    # Analysis of the findings
    print("\nClinical Analysis:")
    if map_pressure < 65:
        print("The MAP of {:.2f} mmHg is below the target of 65 mmHg, indicating potential inadequate organ perfusion.".format(map_pressure))
    else:
        print("The MAP of {:.2f} mmHg is above the critical threshold of 65 mmHg, but the patient's overall clinical picture must be considered.".format(map_pressure))
    
    print("The patient presents with hypotension (90/60), tachycardia (HR 100), tachypnea (RR 40), and dehydration.")
    print("This clinical picture is highly suggestive of shock, likely septic shock given the necrotic tissue.")
    print("The immediate priorities in managing shock are restoring circulating volume and treating the underlying cause.")
    print("A. Intravenous fluids are essential to correct dehydration and improve blood pressure.")
    print("B. Intravenous medications (antibiotics) are critical as PO/topical medications have failed and there is evidence of a severe systemic process.")
    print("C. Surgical debridement is necessary for source control but follows initial resuscitation.")
    print("E. High-flow O2 is not the first priority as SpO2 is 98%.")
    print("Therefore, the most appropriate initial combination of treatments is A and B.")

calculate_map()