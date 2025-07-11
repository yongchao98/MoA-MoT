import textwrap

def present_catalyst_design():
    """
    Presents a hypothetical optimal catalyst system for dual-function
    polyolefin polymerization and hydrogenolysis.
    """

    # --- Catalyst System Definition ---
    catalyst_system = {
        "metal": {
            "name": "Zirconium (Zr)",
            "rationale": [
                "Exhibits a superb balance of high activity for both olefin insertion (polymerization) and sigma-bond metathesis pathways like β-alkyl elimination, which is a key step in polymer breakdown.",
                "Zirconium-based metallocenes and constrained-geometry catalysts are well-established, high-performance polymerization catalysts."
            ]
        },
        "ligand": {
            "name": "Constrained-Geometry Ligand ({η⁵-C₅Me₄(SiMe₂NtBu)}²⁻)",
            "rationale": [
                "Its unique geometry creates a highly electrophilic and sterically accessible metal center, leading to extremely high catalytic activity.",
                "The covalent 'tether' between the cyclopentadienyl and amido groups provides exceptional thermal stability, crucial for the higher temperatures required for hydrogenolysis.",
                "This open structure is key to enabling β-alkyl elimination, the mechanism needed to 'unzip' the polymer chain."
            ]
        },
        "support": {
            "name": "Mesoporous Silica (SBA-15)",
            "rationale": [
                "Provides a very high surface area (>600 m²/g), which allows for the immobilization of the catalyst as isolated single sites.",
                "Single-site isolation prevents bimolecular decomposition pathways, increasing catalyst lifetime and maintaining consistent activity.",
                "Offers excellent mechanical and thermal stability for use in industrial reactors."
            ]
        }
    }

    # --- Introduction ---
    print("### Design for a Dual-Function Single-Site Catalyst ###")
    print("-" * 70)
    print("This is a complex research problem. The following presents a plausible, high-performance catalyst system based on established chemical principles.")
    print("-" * 70)

    # --- Detailed Component Breakdown ---
    for component, details in catalyst_system.items():
        print(f"\n[Component: {component.capitalize()}]")
        print(f"  Suggested Choice: {details['name']}")
        print("  Rationale:")
        # Use textwrap to format the rationale nicely
        for line in details['rationale']:
            wrapped_text = textwrap.fill(line, width=65, initial_indent='    - ', subsequent_indent='      ')
            print(wrapped_text)

    # --- Operating Principle ---
    print("\n[Operating Principle]")
    polymerization_desc = "In the absence of H₂ and presence of olefin monomers (e.g., ethylene), the catalyst rapidly forms high molecular weight polyolefins."
    hydrogenolysis_desc = "At elevated temperature and high H₂ pressure, the catalyst's function reverses. It facilitates tandem β-alkyl elimination and hydrogenation to efficiently cleave the polymer backbone into short-chain alkanes."
    print(textwrap.fill(f"  Polymerization Mode: {polymerization_desc}", width=70, initial_indent='  ', subsequent_indent='  '))
    print(textwrap.fill(f"  Hydrogenolysis Mode: {hydrogenolysis_desc}", width=70, initial_indent='  ', subsequent_indent='  '))
    
    # --- Final "Equation" Representation ---
    print("-" * 70)
    print("\nConceptual 'Equation' for Catalyst Assembly:")
    print("The following represents the components that form the final catalyst. As requested, each component of the 'equation' is explicitly named below.")

    metal_choice = f"1 * Metal Center ({catalyst_system['metal']['name']})"
    ligand_choice = f"1 * Ligand ({catalyst_system['ligand']['name']})"
    support_choice = f"1 * Support ({catalyst_system['support']['name']})"

    print(f"\n{support_choice}")
    print("+")
    print(f"{ligand_choice}")
    print("+")
    print(f"{metal_choice}")
    print("\n--->  [Supported Single-Site Dual-Function Catalyst]")
    print("-" * 70)

# Execute the function to print the analysis
present_catalyst_design()
