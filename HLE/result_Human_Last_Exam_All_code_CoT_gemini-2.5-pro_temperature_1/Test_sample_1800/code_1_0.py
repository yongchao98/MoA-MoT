import textwrap

def get_ni_ce_ratio_info():
    """
    Provides information on the ideal Ni/Ce ratio for Ni-Ceria catalysts.
    """
    explanation = """
    The ideal Ni/Ce atomic ratio in Ni-Ceria nanoparticles for maximizing catalytic performance in the Water Gas Shift (WGS) and Water Splitting (WS) reactions is not a single fixed value, as it depends on synthesis methods and reaction conditions. However, there is a well-established principle based on balancing two key factors:

    1.  **Dispersion of Nickel:** A lower Ni content promotes high dispersion of Ni particles on the Ceria support. This maximizes the crucial interface between Ni and Ceria, which is where the catalytic action occurs.

    2.  **Number of Active Sites:** A higher Ni content provides more active sites for the reaction. However, excessive nickel can lead to the agglomeration (sintering) of Ni particles, which drastically reduces the catalyst's surface area and long-term stability.

    Based on extensive research in catalysis, a consensus has emerged. An optimal balance is typically struck with a relatively low but sufficient amount of nickel.
    """

    conclusion = """
    A highly effective composition that balances high nickel dispersion with a sufficient number of active sites is frequently found to be a Ni/Ce atomic ratio of 1 to 4.
    """

    final_ratio = 0.25

    print("--- Analysis of Ideal Ni/Ce Ratio in Ni-Ceria Catalysts ---")
    print(textwrap.dedent(explanation))
    print(textwrap.dedent(conclusion))
    print("Therefore, the final recommended equation is:")
    print("Ideal Ni / Ce Ratio = 1 / 4")
    print(f"Which corresponds to a numerical value of: {final_ratio}")


if __name__ == "__main__":
    get_ni_ce_ratio_info()