def analyze_clinical_data():
    """
    This function extracts and displays the numerical data from the described
    clinical scenario and formats it into a summary equation.
    """
    # Numerical data extracted from the text
    glycated_hemoglobin = 7.5  # in percent
    time_to_hospital = 28  # in hours
    snb_angle = 88  # in degrees

    # Individual teeth lost
    lost_lateral_incisor = 1
    lost_canine = 1
    lost_first_premolar = 1

    # Calculating total teeth lost for clarity
    total_teeth_lost = lost_lateral_incisor + lost_canine + lost_first_premolar

    print("Extracted Numerical Data from the Scenario:")
    print(f"Glycated Hemoglobin Level: {glycated_hemoglobin}%")
    print(f"Time After Accident: {time_to_hospital} hours")
    print(f"SNB Angle: {snb_angle} degrees")
    print(f"Total Number of Teeth Lost: {total_teeth_lost}")
    print("-" * 30)

    # As requested, printing all numbers in a final equation format
    print("Final Equation with All Numerical Values:")
    # The equation includes each distinct numerical value mentioned in the prompt.
    print(f"{glycated_hemoglobin} + {time_to_hospital} + {snb_angle} + {lost_lateral_incisor} + {lost_canine} + {lost_first_premolar}")


if __name__ == "__main__":
    analyze_clinical_data()