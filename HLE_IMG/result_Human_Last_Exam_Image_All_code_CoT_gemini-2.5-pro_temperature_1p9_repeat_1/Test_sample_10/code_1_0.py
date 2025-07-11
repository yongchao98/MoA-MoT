# Data for the three molecules
molecules = {
    'Corannulene': {'C': 20, 'H': 10, 'barrier': 10},
    'Diindenocorannulene': {'C': 36, 'H': 12, 'barrier': 49},
    'Hemifullerene C42H18': {'C': 42, 'H': 18, 'barrier': None}
}

def calculate_bonds(C, H):
    """Calculates the number of C-C bonds in a PAH."""
    return (3 * C - H) / 2

# Calculate the number of C-C bonds for each molecule
for name, data in molecules.items():
    data['bonds'] = calculate_bonds(data['C'], data['H'])

# Use the data for the first two molecules to build a linear model: barrier = m * bonds + c
m1_data = molecules['Corannulene']
m2_data = molecules['Diindenocorannulene']

E1 = m1_data['bonds']
B1 = m1_data['barrier']
E2 = m2_data['bonds']
B2 = m2_data['barrier']

# Calculate slope (m)
m = (B2 - B1) / (E2 - E1)

# Calculate intercept (c)
c = B1 - m * E1

print(f"--- Model Development ---")
print(f"Molecule 1 (Corannulene): C={m1_data['C']}, H={m1_data['H']} -> Bonds={E1}, Barrier={B1} kcal/mol")
print(f"Molecule 2 (Diindenocorannulene): C={m2_data['C']}, H={m2_data['H']} -> Bonds={E2}, Barrier={B2} kcal/mol")
print("\nAssuming a linear relationship: Barrier = m * (Number of Bonds) + c")
print(f"Slope (m) = ({B2} - {B1}) / ({E2} - {E1}) = {m:.4f}")
print(f"Intercept (c) = {B1} - {m:.4f} * {E1} = {c:.4f}")
print(f"Final Model: Barrier = {m:.4f} * (Number of Bonds) + {c:.4f}\n")


# Predict the barrier for the third molecule
m3_data = molecules['Hemifullerene C42H18']
E3 = m3_data['bonds']

# Calculate the predicted barrier
B3_predicted = m * E3 + c
B3_rounded = round(B3_predicted)

print(f"--- Prediction for the Third Molecule ---")
print(f"Molecule 3 (Hemifullerene): C={m3_data['C']}, H={m3_data['H']} -> Bonds={E3}")
print("Applying the model to predict the inversion barrier:")
print(f"Predicted Barrier = {m:.4f} * {E3} + {c:.4f} = {B3_predicted:.4f} kcal/mol")
print(f"\nRounding to the nearest integer, the predicted inversion barrier is {B3_rounded} kcal/mol.")
print("\nThe equation with the final numbers is:")
print(f"{B3_rounded} = ({m:.4f} * {E3}) + ({c:.4f})")
