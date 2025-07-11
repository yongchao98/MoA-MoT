def solve_insect_wings():
    """
    Analyzes the provided insect wings to determine the trophic level of the family each belongs to.
    """
    
    # Analysis of each wing
    wing_a = {
        "ID": "A",
        "Order": "Hymenoptera",
        "Family": "Ichneumonidae (Ichneumon Wasp)",
        "Venation": "Features a prominent pterostigma and a characteristic areolet cell. The venation is reduced, which is typical for many parasitic wasps.",
        "Trophic Level": "Parasitoid"
    }

    wing_b = {
        "ID": "B",
        "Order": "Neuroptera",
        "Family": "Chrysopidae or Hemerobiidae (Lacewing)",
        "Venation": "Complex, net-like venation with extensive branching of longitudinal veins and numerous crossveins.",
        "Trophic Level": "Predator"
    }
    
    wing_c = {
        "ID": "C",
        "Order": "Hemiptera",
        "Family": "Pentatomidae (Stink Bug)",
        "Venation": "This is a hemelytron, a wing that is hardened at the base (corium) and membranous at the tip. Veins are visible in the membrane, forming a few large cells.",
        "Trophic Level": "Herbivore"
    }

    # Print the detailed description
    print("Detailed Wing Analysis:")
    for wing in [wing_a, wing_b, wing_c]:
        print(f"--- Wing {wing['ID']} ---")
        print(f"Order: {wing['Order']}")
        print(f"Probable Family: {wing['Family']}")
        print(f"Venation Description: {wing['Venation']}")
        print(f"Trophic Level: {wing['Trophic Level']}")
        print("-" * 20)

    # Print the final combined answer
    print("\nFinal Conclusion:")
    print(f"A: {wing_a['Trophic Level']}, B: {wing_b['Trophic Level']}, C: {wing_c['Trophic Level']}")

solve_insect_wings()