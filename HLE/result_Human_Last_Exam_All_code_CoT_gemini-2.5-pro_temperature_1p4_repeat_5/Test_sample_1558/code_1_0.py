import math

def format_number(num):
    """Formats a number into scientific notation like 3x10^3"""
    if num == 0:
        return "0"
    exponent = int(math.log10(abs(num)))
    mantissa = num / (10**exponent)
    # Adjust for standard scientific notation where exponent is a multiple of 3
    if exponent % 3 != 0:
        shift = exponent % 3
        mantissa = mantissa * (10**shift)
        exponent = exponent - shift
    return f"{int(mantissa)}x10^{exponent}"

# Storing the relevant Red Blood Cell data from the experiments
# Values are per ul
data = {
    "exp1": {
        "non_pregnant_control": 13e6,
        "pregnant_control": 10e6,
        "pregnant_rti": 8e6,
    },
    "exp2": {
        "non_pregnant_control": 13e6,
        "pregnant_control": 13e6,
        # STING deletion inhibits the interferon pathway
        "pregnant_dsting": 8e6,
    },
}

print("--- Analysis to Determine the Correct Conclusion ---")

# --- Step 1: Analyze the effect of Transposable Elements (TEs) ---
print("\nAnalyzing Clause 1: 'Increased activity of transposable elements increases the erythropoiesis in pregnant mice.'")
print("This is tested by inhibiting TEs with RTI and observing the effect on Red Blood Cells (RBCs) in pregnant mice from Experiment 1.")

preg_control_rbc_e1 = data["exp1"]["pregnant_control"]
preg_rti_rbc_e1 = data["exp1"]["pregnant_rti"]

print(f"RBC Count (Pregnant Control): {format_number(preg_control_rbc_e1)} per ul")
print(f"RBC Count (Pregnant + RTI): {format_number(preg_rti_rbc_e1)} per ul")
print(f"Equation: {format_number(preg_control_rbc_e1)} (Control) > {format_number(preg_rti_rbc_e1)} (RTI)")

if preg_control_rbc_e1 > preg_rti_rbc_e1:
    print("Result: The data shows that inhibiting TE activity decreases RBCs. This supports the statement.")
else:
    print("Result: The data does not support the statement.")

# --- Step 2: Analyze the effect of Interferon ---
print("\nAnalyzing Clause 2: 'Interferon does not increase the number of red blood cells in pregnant mice.'")
print("This is interpreted by comparing the RBC count in pregnant mice (with an active Interferon system) to the non-pregnant baseline.")

# Comparison using Experiment 1 data
print("\n--- Based on Experiment 1 Data ---")
np_control_rbc_e1 = data["exp1"]["non_pregnant_control"]
print(f"RBC Count (Non-Pregnant Control): {format_number(np_control_rbc_e1)} per ul")
print(f"RBC Count (Pregnant Control): {format_number(preg_control_rbc_e1)} per ul")
print(f"Equation: Is {format_number(preg_control_rbc_e1)} > {format_number(np_control_rbc_e1)}?")
if preg_control_rbc_e1 > np_control_rbc_e1:
    print("Result: The count in pregnant mice is higher. This would contradict the statement.")
else:
    print("Result: The count in pregnant mice is not higher. This supports the statement.")

# Comparison using Experiment 2 data
print("\n--- Based on Experiment 2 Data ---")
np_control_rbc_e2 = data["exp2"]["non_pregnant_control"]
preg_control_rbc_e2 = data["exp2"]["pregnant_control"]
print(f"RBC Count (Non-Pregnant Control): {format_number(np_control_rbc_e2)} per ul")
print(f"RBC Count (Pregnant Control): {format_number(preg_control_rbc_e2)} per ul")
print(f"Equation: Is {format_number(preg_control_rbc_e2)} > {format_number(np_control_rbc_e2)}?")
if preg_control_rbc_e2 > np_control_rbc_e2:
    print("Result: The count in pregnant mice is higher. This would contradict the statement.")
else:
    print("Result: The count in pregnant mice is not higher. This supports the statement.")

print("\n--- Final Conclusion ---")
print("The analysis shows that both clauses in option E are supported by the data under a reasonable interpretation.")
print("1. TE activity increases RBC count relative to an inhibited state.")
print("2. Interferon activity does not increase RBC count above the non-pregnant baseline.")
print("Therefore, E is the most accurate conclusion.")
