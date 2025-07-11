import math

# Store the experimental data
data = {
    "Control": 500,
    "Zma1 + standard buffer + 5 mM MgCl2": 700,
    "Zma1 + standard buffer + 5 mM CaCl2": 500,
    "Zma1 + standard buffer + 5 mM CuCl2": 400,
    "Zma1 + standard buffer + 5 mM Al1": 1000,
    "Zma1 + standard buffer + 5mM Al2": 150,
    "Zma1 + standard buffer + 5mM Al1 + 5mM Al2": 150,
    "Zma1 + standard buffer + 100 mM XAG1": 10,
    "Zma1 + standard buffer + 100 mM XAG1 + 500 mM of molecule A": 450,
    "Zma1 + standard buffer + 100 mM Rga1": 10,
    "Zma1 + standard buffer + 100 mM Rga1 + 500 mM of molecule A": 10,
}

kcat_control = data["Control"]

print("--- Analysis of Enzyme Kinetics Data for Zma1 ---")
print(f"1. The baseline (control) activity of Zma1 is kcat = {kcat_control}/second.\n")

print("--- Analyzing the function of Al1 ---")
kcat_al1 = data["Zma1 + standard buffer + 5 mM Al1"]
print(f"Adding Al1 increases kcat from {kcat_control} to {kcat_al1}/second.")
change_factor_al1 = kcat_al1 / kcat_control
print(f"This is a {change_factor_al1}-fold increase in activity.")
print("Conclusion: Al1 is an allosteric activator.\n")

print("--- Analyzing the function of Rga1 ---")
kcat_rga1 = data["Zma1 + standard buffer + 100 mM Rga1"]
kcat_rga1_high_A = data["Zma1 + standard buffer + 100 mM Rga1 + 500 mM of molecule A"]
print(f"Adding Rga1 decreases kcat from {kcat_control} to {kcat_rga1}/second.")
print(f"When excess substrate (Molecule A) is added, the kcat remains at {kcat_rga1_high_A}/second.")
print("Since high substrate concentration does not reverse the inhibition, Rga1 is not a competitive inhibitor.")
print("Conclusion: Rga1 is a non-competitive or irreversible inhibitor.\n")

print("--- Final Evaluation based on Analysis ---")
print("Based on the data:")
print("- Al1 functions as an allosteric activator.")
print("- Rga1 functions as a non-competitive or irreversible inhibitor.")
print("- Additional analysis shows Al1 and Al2 compete for the same site because with both present, the kcat is {data['Zma1 + standard buffer + 5mM Al1 + 5mM Al2']}/s, identical to Al2 alone ({data['Zma1 + standard buffer + 5mM Al2']}/s).")

print("\nEvaluating the choices, option C aligns best with all the evidence:")
print("C. Al1 and Al2 function as allosteric modulators for Zma1. Al1 and Al2 bind to the same site on Zma1. Rga1 is an irreversible inhibitor.")
<<<C>>>