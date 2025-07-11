def calculate_heritability(va, vd, vi, ve):
    """
    Calculates broad-sense (H2) and narrow-sense (h2) heritability
    from variance components.
    va: Additive genetic variance
    vd: Dominance genetic variance
    vi: Epistatic genetic variance
    ve: Environmental variance
    """
    # Calculate total genetic and phenotypic variances
    vg = va + vd + vi
    vp = vg + ve

    # Avoid division by zero
    if vp == 0:
        return 0, 0

    # Calculate heritabilities
    h2_broad = vg / vp
    h2_narrow = va / vp
    return h2_broad, h2_narrow, vg, vp

# --- Scenario 1: Rabbits (Purely Additive Genetics) ---
# Given H^2 = 0.75. We can set VA=75 and VE=25 to model this.
# For rabbits, VD and VI are 0.
va_rabbit = 75
vd_rabbit = 0
vi_rabbit = 0
ve_rabbit = 25

H2_rab, h2_rab, Vg_rab, Vp_rab = calculate_heritability(va_rabbit, vd_rabbit, vi_rabbit, ve_rabbit)

print("--- Rabbit Scenario Analysis ---")
print(f"Given: VA={va_rabbit}, VD={vd_rabbit}, VI={vi_rabbit}, VE={ve_rabbit}")
print(f"Equation for H^2: VG / VP = ({va_rabbit} + {vd_rabbit} + {vi_rabbit}) / {Vp_rab} = {Vg_rab} / {Vp_rab} = {H2_rab:.2f}")
print(f"Equation for h^2: VA / VP = {va_rabbit} / {Vp_rab} = {h2_rab:.2f}")
print("Result: When genetic variance is purely additive, H^2 and h^2 are equal.\n")


# --- Scenario 2: Testing the Answer Choices ---
print("--- Evaluating Answer Choices ---")

# Test Choice E: Does Dominance Variance (VD) have no impact on h^2?
# Let's introduce VD=20 to the rabbit model.
vd_species = 20
H2_E, h2_E, Vg_E, Vp_E = calculate_heritability(va_rabbit, vd_rabbit + vd_species, vi_rabbit, ve_rabbit)
print("[Analysis of E] Effect of adding Dominance Variance:")
print(f"With VD={vd_species}, the new h^2 equation is: VA / VP = {va_rabbit} / {Vp_E} = {h2_E:.2f}")
print(f"Conclusion: The statement in E is FALSE. Adding VD changed h^2 from {h2_rab:.2f} to {h2_E:.2f}.\n")


# Test Choice C: Does Epistatic Interaction (VI) cause measures to differ?
# Let's introduce VI=20 to the rabbit model.
vi_species = 20
H2_C, h2_C, Vg_C, Vp_C = calculate_heritability(va_rabbit, vd_rabbit, vi_rabbit + vi_species, ve_rabbit)
print("[Analysis of C] Effect of adding Epistatic Variance:")
print(f"With VI={vi_species}, the new H^2 equation is: VG / VP = {Vg_C} / {Vp_C} = {H2_C:.2f}")
print(f"With VI={vi_species}, the new h^2 equation is: VA / VP = {va_rabbit} / {Vp_C} = {h2_C:.2f}")
print(f"Conclusion: The statement in C is TRUE. The presence of epistasis caused H^2 ({H2_C:.2f}) and h^2 ({h2_C:.2f}) to differ from each other, and it caused the value of h^2 to change (from {h2_rab:.2f} to {h2_C:.2f}).")
