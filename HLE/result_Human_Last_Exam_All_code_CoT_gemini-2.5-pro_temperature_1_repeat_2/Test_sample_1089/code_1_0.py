def main():
    # The problem asks to identify the main source of damage on aeroengine turbine blades
    # that is repaired by manual TIG welding build-up.

    # The repair method is key: "TIG welding (GTAW) build-up of layers of filler material".
    # This is an additive process, meaning it's used to replace material that has been lost.

    # Let's consider the answer choices:
    # A. Stress Corrosion Cracking - A crack, not primarily a material loss issue.
    # B. Foreign Object Damage - Involves material loss (nicks, gouges), a possible answer.
    # C. Blade Tip Rub and Wear - Involves material being worn away from the blade tip. This is a very common issue.
    # D. Creep Deformation - Blade stretching/lengthening; adding material is not the fix.
    # E. Fatigue Cracking - A crack, similar to A.
    # F. High-Temperature Oxidation and Corrosion - Involves surface material loss, another possibility.

    # Comparing the material loss options (B, C, F):
    # Blade Tip Rub and Wear (C) is a systematic and frequent issue in engine operation.
    # Restoring the precise blade tip clearance is critical for engine efficiency.
    # The process of building up layers of weld material on the tip and then re-machining it to the
    # original specifications is the classic application of this repair technique in MRO.
    # It directly addresses the loss of "geometrical integrity" at the blade tip.

    # While FOD and corrosion can also be repaired with welding, the most common and direct
    # application of a "build-up of layers" to restore geometry is for tip wear.

    answer = 'C'
    explanation = f"The main source of damage addressed by manual TIG welding build-up is Blade Tip Rub and Wear. This is because the repair method is additive, designed to replace lost material and restore the blade's original geometry. Blade tip wear is a common issue where material is ground off the tip, and TIG welding is the standard procedure to build it back up to maintain engine efficiency."

    print(f"The correct choice is {answer}.")
    print("\nExplanation:")
    print(explanation)

if __name__ == "__main__":
    main()