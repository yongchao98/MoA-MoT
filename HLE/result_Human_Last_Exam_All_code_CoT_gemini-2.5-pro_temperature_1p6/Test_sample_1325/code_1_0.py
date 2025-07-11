def calculate_heritability(Va, Vd, Vi, Ve):
    """Calculates and prints heritability values based on variance components."""
    Vg = Va + Vd + Vi
    Vp = Vg + Ve

    # Prevent division by zero
    if Vp == 0:
        print("Phenotypic variance is zero, cannot calculate heritability.")
        return

    H2 = Vg / Vp
    h2 = Va / Vp
    
    print(f"Variance Components: Va={Va}, Vd={Vd}, Vi={Vi}, Ve={Ve}")
    print(f"Total Genotypic Variance (Vg) = Va + Vd + Vi = {Va} + {Vd} + {Vi} = {Vg}")
    print(f"Phenotypic Variance (Vp) = Vg + Ve = {Vg} + {Ve} = {Vp}")
    print("-" * 20)
    print("Broad-Sense Heritability (H^2) = Vg / Vp")
    print(f"H^2 = {Vg} / {Vp} = {H2:.2f}")
    print("\nNarrow-Sense Heritability (h^2) = Va / Vp")
    print(f"h^2 = {Va} / {Vp} = {h2:.2f}")
    
    if H2 == h2:
        print("\nResult: H^2 and h^2 are equal.")
    else:
        print("\nResult: H^2 and h^2 are different.")


# --- Scenario 1: Rabbits with purely additive genetics ---
# The problem states H^2 = 0.75. If genetics are purely additive, H^2 = h^2 = Va / (Va + Ve).
# So, 0.75 = Va / (Va + Ve) => 0.75 * Ve = 0.25 * Va => Ve = Va / 3.
# Let's choose Va = 30, so Ve = 10.
print("--- Scenario 1: Rabbit Population (Purely Additive Genetics) ---")
# Here Vd and Vi are 0 as the model is purely additive.
calculate_heritability(Va=30, Vd=0, Vi=0, Ve=10)

print("\n" + "="*60 + "\n")

# --- Scenario 2: Another species with epistatic interactions (Choice C) ---
# We introduce epistatic variance (Vi > 0) to see the effect.
print("--- Scenario 2: Another Species (with Epistatic Interactions) ---")
# Let's add Vi = 5, keeping other components the same for comparison.
calculate_heritability(Va=30, Vd=0, Vi=5, Ve=10)
