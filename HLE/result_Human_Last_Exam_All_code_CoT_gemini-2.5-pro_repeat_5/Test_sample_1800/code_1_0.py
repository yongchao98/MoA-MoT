import textwrap

def get_ideal_ni_ce_ratio_info():
    """
    This script explains the ideal Ni/Ce ratio for Ni-Ceria catalysts
    in Water Gas Shift (WGS) and Water Splitting (WS) reactions based on
    scientific consensus.
    """

    print("--- Analysis of Ideal Ni/Ce Ratio for Ni-Ceria Catalysts ---")

    # Principle of catalysis
    principle = """
    The catalytic performance of Ni-Ceria nanoparticles for the Water Gas Shift (WGS) and water splitting reactions is primarily determined by the active sites at the interface between Ni and the CeO2 support. The ideal ratio aims to maximize these active sites.
    """
    print("\n[Principle]")
    print(textwrap.fill(principle, width=80))

    # Explanation of the trade-off
    trade_off = """
    There is a critical trade-off in determining the ratio. A low Ni content (<0.1 molar ratio) ensures high dispersion of Ni particles but may not provide enough total active sites. A high Ni content (>0.3 molar ratio) often leads to the sintering of Ni into larger, less effective particles, which reduces the active interface area.
    """
    print("\n[The Trade-off]")
    print(textwrap.fill(trade_off, width=80))

    # Conclusion on the optimal ratio
    conclusion = """
    Based on numerous studies, the optimal performance for both WGS and water splitting is typically achieved with a Ni/Ce molar ratio in a specific range. This range balances the need for sufficient Ni sites with the requirement of high dispersion and a strong Ni-support interaction.
    """
    print("\n[Conclusion]")
    print(textwrap.fill(conclusion, width=80))

    # Final Answer
    # The 'equation' here is the statement of the final optimal range.
    print("\n--- Final Recommended Ratio ---")
    print("The ideal Ni/Ce molar ratio is generally found to be in the range of:")
    lower_bound = 0.1
    upper_bound = 0.3
    print(f"{lower_bound} to {upper_bound}")
    print("\nThis means for every 1 atom of Nickel (Ni), there should ideally be between 3 to 10 atoms of Cerium (Ce).")


if __name__ == "__main__":
    get_ideal_ni_ce_ratio_info()
