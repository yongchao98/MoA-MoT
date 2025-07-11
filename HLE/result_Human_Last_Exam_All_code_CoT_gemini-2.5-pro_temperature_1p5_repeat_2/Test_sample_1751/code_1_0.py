def explain_lipid_packing():
    """
    Explains which lipid, C16-dihydroceramide or C16-ceramide,
    has a lower surface area in a compressed monolayer and why.
    """
    
    lipid_1 = "C16-dihydroceramide"
    lipid_2 = "C16-ceramide"

    print("Analyzing the molecular structures and their effect on packing in a monolayer:")
    print("-" * 70)

    # Step 1: Analyze C16-dihydroceramide
    print(f"1. {lipid_1} (d18:0/16:0):")
    print("   - This lipid is composed of two fully saturated hydrocarbon chains.")
    print("   - Saturated chains are straight and flexible, allowing them to pack together very tightly and orderly.")
    print("   - This tight packing maximizes van der Waals interactions.")
    print()

    # Step 2: Analyze C16-ceramide
    print(f"2. {lipid_2} (d18:1/16:0):")
    print("   - This lipid contains a trans double bond in its sphingoid base (d18:1).")
    print("   - This double bond introduces a rigid point in the chain, disrupting the perfect linear alignment possible with saturated chains.")
    print("   - This disruption hinders the ability of the molecules to pack as tightly as their fully saturated counterparts.")
    print()

    # Step 3: Relate packing to surface area
    print("3. Packing and Surface Area:")
    print("   - When a lipid monolayer is compressed, the molecules are forced into their most compact arrangement.")
    print("   - The surface area per molecule is determined by how tightly the molecules can pack.")
    print("   - Tighter packing results in a smaller surface area per molecule.")
    print()
    
    # Step 4: Conclusion
    print("Conclusion:")
    print(f"Because the fully saturated chains of {lipid_1} allow for tighter, more ordered packing than the chains of {lipid_2} (which are disrupted by a double bond), {lipid_1} will occupy less space upon compression.")
    print("\nTherefore, the lipid with the lower surface area when compressed is:")
    print(f"--> {lipid_1}")


if __name__ == "__main__":
    explain_lipid_packing()