import sys

def analyze_mfa_requirements():
    """
    Analyzes and prints the requirements for a 13C metabolic flux analysis at steady state.
    """

    # List of information items. Each tuple contains:
    # (Description, Is_Required_Boolean, Explanation_If_Required, Explanation_If_Not_Required)
    information_list = [
        (
            "Metabolic reaction stoichiometry",
            True,
            "This forms the fundamental model of the metabolic network. It defines all possible pathways for carbon flow.",
            ""
        ),
        (
            "Maximum cell density of the organism in a bioreactor",
            False,
            "",
            "This is a measure of growth performance, not a direct input for calculating the relative flux distribution at steady state, which is the core of MFA."
        ),
        (
            "Biomass composition of the organism (i.e. fraction of protein, lipids, and carbohydrates)",
            True,
            "This is essential to quantify the metabolic demands for creating new cell components. It provides critical constraints on the fluxes drained from central metabolism.",
            ""
        ),
        (
            "Enzyme kinetics",
            False,
            "",
            "This information is required for dynamic or kinetic modeling, not for steady-state MFA, which calculates the net flux without considering the underlying enzyme mechanics."
        ),
        (
            "Regulatory networks",
            False,
            "",
            "Steady-state MFA measures the outcome of regulation (the flux distribution), but the regulatory network itself is not an input for the calculation."
        ),
        (
            "Isotope labeling patterns of metabolic intermediates",
            True,
            "This is the defining experimental data for 13C-MFA. The labeling patterns provide powerful constraints that allow for the precise calculation of intracellular fluxes.",
            ""
        )
    ]

    required_items_indices = []
    
    # Use ANSI escape codes for color if the terminal supports it
    is_tty = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()
    green_color = "\033[92m" if is_tty else ""
    red_color = "\033[91m" if is_tty else ""
    end_color = "\033[0m" if is_tty else ""

    print("Analyzing the requirements for 13C Metabolic Flux Analysis (MFA) at steady state:\n")

    for i, (item, is_required, req_exp, not_req_exp) in enumerate(information_list, 1):
        if is_required:
            required_items_indices.append(str(i))
            print(f"{i}. {item}: {green_color}Required{end_color}")
            print(f"   - Why: {req_exp}")
        else:
            print(f"{i}. {item}: {red_color}Not Required{end_color}")
            print(f"   - Why: {not_req_exp}")
        print("-" * 70)

    required_count = len(required_items_indices)
    
    print("\n" + "="*28 + " SUMMARY " + "="*29)
    print(f"The analysis shows that {required_count} of the {len(information_list)} listed items are required.")
    print(f"The numbers of the required items are: {', '.join(required_items_indices)}.")
    print("=" * 65)

if __name__ == '__main__':
    analyze_mfa_requirements()
