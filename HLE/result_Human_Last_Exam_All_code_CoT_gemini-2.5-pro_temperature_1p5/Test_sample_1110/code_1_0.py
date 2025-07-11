import sys

def calculate_reaction_outcome(buLi_equivalents):
    """
    Models the formation of a boronic acid and a diarylborinic acid byproduct.

    Args:
        buLi_equivalents (float): The equivalents of n-BuLi used relative to the starting material.
    """
    # Let's assume we start with 100 mmol of the aryl iodide.
    aryl_iodide_moles = 100.0

    # Calculate the moles of n-BuLi added.
    nBuLi_moles = aryl_iodide_moles * buLi_equivalents

    # The problem is caused by excess nucleophile attacking the desired product.
    # Let's model this: the amount of byproduct formed is proportional to the excess n-BuLi.
    # Let's set a simple proportionality constant.
    byproduct_formation_factor = 2.0

    print(f"--- Running simulation with {buLi_equivalents} eq. of n-BuLi ---")

    if buLi_equivalents > 1.0:
        # Excess n-BuLi is present after all aryl iodide is converted.
        excess_nBuLi_moles = nBuLi_moles - aryl_iodide_moles
        
        # Calculate the amount of byproduct formed.
        # This byproduct consumes the desired product (boronic ester).
        byproduct_moles = byproduct_formation_factor * excess_nBuLi_moles
        
        # The amount of desired product is what's left over.
        # It's initially equal to the starting material, then some is consumed.
        desired_product_moles = aryl_iodide_moles - byproduct_moles
        
        # Ensure moles are not negative
        desired_product_moles = max(0, desired_product_moles)
        
        print(f"Excess n-BuLi: {excess_nBuLi_moles:.1f} mmol")
        print("An excess of nucleophile leads to a side reaction, forming a byproduct.")
        print("Final Equation:")
        print(f"Desired Product Moles = {desired_product_moles:.1f}")
        print(f"Byproduct Moles = {byproduct_moles:.1f}")
        print("Result: Two boron-containing products are formed, leading to 2 NMR signals.\n")
        
    else:
        # If n-BuLi is the limiting reagent or stoichiometric, no excess exists.
        # The reaction yield is limited by the amount of n-BuLi.
        desired_product_moles = nBuLi_moles
        byproduct_moles = 0.0
        
        print("No excess n-BuLi is present.")
        print("The side reaction forming the byproduct is suppressed.")
        print("Final Equation:")
        print(f"Desired Product Moles = {desired_product_moles:.1f}")
        print(f"Byproduct Moles = {byproduct_moles:.1f}")
        print("Result: Only one boron-containing product is formed, leading to 1 NMR signal.\n")


# Case from the problem description
calculate_reaction_outcome(1.05)

# Proposed solution
calculate_reaction_outcome(1.00)
