import sys

# Step 1: Define the compound, target, and known binding affinity value.
# The compound 7-dimethylphosphoryl-3-[...]-carbonitrile is also known as SY-1365.
# The target is Cyclin-dependent kinase 7 (CDK7).
# Data is sourced from scientific literature (e.g., Mol Cancer Ther. 2017;16(12 Suppl):Abstract nr PR06).
compound_name = "SY-1365"
target_protein = "Cyclin-dependent kinase 7 (CDK7)"
activity_value = 1.4  # The IC50 value
activity_unit = "nM"  # The unit for the IC50 value

print(f"Finding the binding affinity of {compound_name} to {target_protein}.")
print(f"Based on published biochemical assays, the IC50 value is {activity_value} {activity_unit}.")
print("-" * 30)

# Step 2: Define the answer choices as numerical ranges in nM.
# Note: 1 uM = 1,000 nM; 1 mM = 1,000,000 nM
answer_choices = {
    "A": {"min": -sys.float_info.max, "max": 0.1, "desc": "< 0.1 nM"},
    "B": {"min": 0.1, "max": 100, "desc": "0.1 - 100 nM"},
    "C": {"min": 100, "max": 100000, "desc": "0.1 - 100 uM (100 - 100,000 nM)"},
    "D": {"min": 100000, "max": 100000000, "desc": "0.1 - 100 mM (100,000 - 100,000,000 nM)"},
    "E": {"min": 100000000, "max": sys.float_info.max, "desc": "> 100 mM (> 100,000,000 nM)"}
}

# Step 3: Determine which category the value falls into.
final_answer = None
for choice, range_info in answer_choices.items():
    # We use a strict less than for the lower bound of choice A, and inclusive for the rest.
    is_in_range = False
    if choice == "A":
        if activity_value < range_info["max"]:
            is_in_range = True
    else:
        if range_info["min"] <= activity_value <= range_info["max"]:
            is_in_range = True

    if is_in_range:
        final_answer = choice
        print(f"Checking range {choice} ({range_info['desc']}):")
        if choice == "A":
             print(f"Is {activity_value} < {range_info['max']}? Yes.")
        else:
             print(f"Is {activity_value} >= {range_info['min']} and {activity_value} <= {range_info['max']}? Yes.")
        print(f"Conclusion: The value {activity_value} {activity_unit} falls into category {choice}.")
        break
    else:
        # Show why it doesn't match for clarity
        print(f"Checking range {choice} ({range_info['desc']}):")
        if choice == "A":
             print(f"Is {activity_value} < {range_info['max']}? No.")
        else:
             print(f"Is {activity_value} >= {range_info['min']} and {activity_value} <= {range_info['max']}? No.")

if not final_answer:
    print("\nThe value does not fall into any of the provided categories.")
else:
    print("\nFinal Answer Category:", final_answer)
