# --- User Data ---
# Positive relationships (Trust Incoming)
p1 = {'total': 7, 'trust': 6, 'distrust': 1}
p2 = {'total': 6, 'trust': 4, 'distrust': 2}
p3 = {'total': 4, 'trust': 2, 'distrust': 2}

# Negative relationships (Distrust Incoming)
n1 = {'total': 6, 'trust': 3, 'distrust': 3}
n2 = {'total': 4, 'trust': 1, 'distrust': 3}

# --- FAT Rule Calculations ---

# Contribution from P1 (Positive Edge)
# Rule 1: 1 / (total_relationships + 1)
p1_score = 1 / (p1['total'] + 1)
p1_eq = f"1/({p1['total']}+1)"

# Contribution from P2 (Positive Edge)
# Rule 1: 1 / (total_relationships + 1)
p2_score = 1 / (p2['total'] + 1)
p2_eq = f"1/({p2['total']}+1)"

# Contribution from P3 (Positive Edge)
# Rule 1: 1 / (total_relationships + 1)
p3_score = 1 / (p3['total'] + 1)
p3_eq = f"1/({p3['total']}+1)"

# Contribution from N1 (Negative Edge with Mixed Ratings)
# Rule 2: -1 / (total+1) * (trust/total)
n1_score = -1 / (n1['total'] + 1) * (n1['trust'] / n1['total'])
n1_eq = f"-1/({n1['total']}+1) * ({n1['trust']}/{n1['total']})"


# Contribution from N2 (Negative Edge with More Distrust, hence 1.5x weight)
# Rule 2 and Rule 3: [-1 / (total+1) * (trust/total)] * 1.5
n2_score = (-1 / (n2['total'] + 1) * (n2['trust'] / n2['total'])) * 1.5
n2_eq = f"(-1/({n2['total']}+1) * ({n2['trust']}/{n2['total']}))*1.5"

# --- Total Score Calculation ---
total_score = p1_score + p2_score + p3_score + n1_score + n2_score

# --- Output Results ---
print("B's Importance Score Calculation (FAT measure):")
print(f"Contribution from P1: {p1_eq} = {p1_score:.4f}")
print(f"Contribution from P2: {p2_eq} = {p2_score:.4f}")
print(f"Contribution from P3: {p3_eq} = {p3_score:.4f}")
print(f"Contribution from N1: {n1_eq} = {n1_score:.4f}")
print(f"Contribution from N2: {n2_eq} = {n2_score:.4f}")
print("\nFinal Equation:")
# Showing the plus/minus signs explicitly for clarity
print(f"Score = {p1_eq} + {p2_eq} + {p3_eq} + ({n1_eq}) + ({n2_eq})")
print(f"\nTotal Score = {p1_score:.4f} + {p2_score:.4f} + {p3_score:.4f} {n1_score:.4f} {n2_score:.4f}")
print(f"Total Score = {total_score:.4f}")

# Compare to options
options = {'A': 0.35, 'B': 0.45, 'C': 0.55, 'D': 0.64}
closest_option = min(options, key=lambda k: abs(options[k] - total_score))
print(f"\nThe calculated score {total_score:.4f} is closest to option {closest_option}: {options[closest_option]}")
print(f"\nFinal Answer: {closest_option}")