# --- Hypothetical Fitness Data (e.g., growth rate in doublings/hour) ---
# Fitness of Parental Line 1 (no gene flow)
fitness_p1 = 0.50
# Fitness of Parental Line 2 (no gene flow)
fitness_p2 = 0.48
# Fitness of the F1 hybrid (direct result of gene flow)
# We'll simulate a case of hybrid vigor where the F1 is slightly fitter
fitness_f1_hybrid = 0.52
# Fitness of the F2 hybrid (after meiosis and recombination)
# Here, we simulate a cost due to the breakdown of co-adapted gene complexes
fitness_f2_hybrid = 0.35

print("This script calculates the selection coefficient (s) as a measure of fitness cost due to gene flow.")
print("A positive 's' indicates a cost, a negative 's' indicates a benefit (hybrid vigor).\n")

# --- Step 1: Calculate the baseline fitness of the 'no gene flow' lines ---
fitness_parent_avg = (fitness_p1 + fitness_p2) / 2
print(f"1. Baseline fitness (average of parents):")
print(f"   W_parent_avg = ({fitness_p1} + {fitness_p2}) / 2 = {fitness_parent_avg:.4f}\n")


# --- Step 2: Calculate the selection coefficient for the F1 hybrid ---
# Relative fitness (w) = W_hybrid / W_parent
# Selection coefficient (s) = 1 - w
relative_fitness_f1 = fitness_f1_hybrid / fitness_parent_avg
s_f1 = 1 - relative_fitness_f1

print(f"2. Cost/benefit for F1 Hybrid (initial gene flow):")
print(f"   s_f1 = 1 - (W_f1_hybrid / W_parent_avg)")
print(f"   s_f1 = 1 - ({fitness_f1_hybrid} / {fitness_parent_avg:.4f})")
print(f"   s_f1 = {s_f1:.4f}")
print("   (A negative value indicates the F1 hybrid is fitter than the parents - hybrid vigor)\n")

# --- Step 3: Calculate the selection coefficient for the F2 hybrid (after meiosis) ---
# This step addresses the "within mating to account for effects of meiosis" part of the answer
relative_fitness_f2 = fitness_f2_hybrid / fitness_parent_avg
s_f2 = 1 - relative_fitness_f2

print(f"3. Cost/benefit for F2 Hybrid (after meiosis):")
print(f"   s_f2 = 1 - (W_f2_hybrid / W_parent_avg)")
print(f"   s_f2 = 1 - ({fitness_f2_hybrid} / {fitness_parent_avg:.4f})")
print(f"   s_f2 = {s_f2:.4f}")
print("   (A positive value indicates the F2 hybrid is less fit, showing a cost of gene flow)")
