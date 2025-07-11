import textwrap

def explain_surface_area():
    """
    Explains why C16-dihydroceramide has a lower surface area than C16-ceramide
    in a compressed monolayer.
    """

    # Define the lipids
    lipid_1 = "C16-dihydroceramide (d18:0/16:0)"
    lipid_2 = "C16-ceramide (d18:1/16:0)"

    # Print explanation
    print("Step 1: Identify the key structural difference.")
    explanation_1 = (f"The primary difference is that {lipid_1} has two fully saturated hydrocarbon chains. "
                     f"In contrast, {lipid_2} has one saturated chain and one unsaturated chain containing a "
                     f"'trans' double bond in its sphingosine base.")
    print(textwrap.fill(explanation_1, width=80))
    print("-" * 80)

    print("Step 2: Relate structure to molecular packing.")
    explanation_2 = ("Saturated hydrocarbon chains are straight and flexible, which allows them to pack together "
                     "very tightly and form ordered structures. The 'trans' double bond in C16-ceramide introduces a "
                     "rigid kink in the chain, which disrupts this tight, orderly packing.")
    print(textwrap.fill(explanation_2, width=80))
    print("-" * 80)
    
    print("Step 3: Connect packing to surface area.")
    explanation_3 = ("When a lipid monolayer is compressed, the total surface area is determined by the minimum area "
                     "each molecule can occupy. Tightly packed molecules, like C16-dihydroceramide, occupy a smaller "
                     "cross-sectional area. The disrupted packing of C16-ceramide forces its molecules to occupy a "
                     "larger area.")
    print(textwrap.fill(explanation_3, width=80))
    print("-" * 80)

    print("Conclusion:")
    conclusion = (f"Therefore, because of its ability to pack more tightly, {lipid_1} will have a lower surface "
                  "area when compressed in a monolayer.")
    print(textwrap.fill(conclusion, width=80))

if __name__ == '__main__':
    explain_surface_area()