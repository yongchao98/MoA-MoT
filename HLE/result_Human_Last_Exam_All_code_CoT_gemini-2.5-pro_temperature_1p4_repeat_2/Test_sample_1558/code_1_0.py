import locale

def analyze_hematopoiesis_data():
    """
    Analyzes data from Experiment 1 to determine the role of transposable elements
    in red blood cell production in pregnant mice.
    """
    # Data from Experiment 1: Red Blood Cells in pregnant mice
    rbc_pregnant_control_val = 10 * 10**6
    rbc_pregnant_rti_val = 8 * 10**6

    # For clear output formatting
    locale.setlocale(locale.LC_ALL, '')

    print("Step 1: State the relevant data from Experiment 1.")
    print(f"Red Blood Cell count in pregnant control mice: {int(rbc_pregnant_control_val):n} per ul.")
    print(f"Red Blood Cell count in pregnant mice treated with RTI (inhibitor of transposable elements): {int(rbc_pregnant_rti_val):n} per ul.")
    print("-" * 50)

    print("Step 2: Calculate the change in red blood cells due to RTI treatment.")
    reduction_val = rbc_pregnant_control_val - rbc_pregnant_rti_val
    
    # Final equation with each number printed
    print("Equation of the change:")
    print(f"{int(rbc_pregnant_control_val):n} (control) - {int(rbc_pregnant_rti_val):n} (RTI treated) = {int(reduction_val):n} (reduction)")
    print("-" * 50)

    print("Step 3: Interpret the results and draw a conclusion.")
    print("The data shows that inhibiting transposable elements reduces the red blood cell count in pregnant mice.")
    print("This implies that the activity of transposable elements normally helps increase or maintain the level of red blood cells (erythropoiesis) in this model.")
    print("\nAnemia is a condition characterized by a low count of red blood cells.")
    print("Since activating transposable elements appears to boost red blood cell production, it is a logical inference that their induction could be explored as a potential treatment for anemia.")

# Run the analysis
analyze_hematopoiesis_data()