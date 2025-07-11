import math

def find_computational_factor():
    """
    This function presents the solution by identifying the computational factor
    from the relevant scientific literature as described in the problem.
    """

    # The problem describes the enthalpy-porosity method for modeling solidification/melting.
    # The two key papers are:
    # 1. Prior simulation-only work: Voller, V. R., & Prakash, C. (1987). "A fixed-grid numerical
    #    modelling methodology for convection-diffusion mushy region phase-change problems."
    #    International Journal of Heat and Mass Transfer, 30(8), 1709-1719.
    # 2. Later benchmarking work: Brent, A. D., Voller, V. R., & Reid, K. J. (1988).
    #    "An enthalpy-porosity technique for modeling convection-diffusion phase change:
    #    Application to the melting of a pure metal." Numerical Heat Transfer, 13(3), 297-318.

    # From the prior (1987) simulation-only paper, the mushy zone constant 'C'
    # (the "computational factor") is stated on page 1711.
    c_prior = 1.6 * 10**6

    # From the later (1988) benchmarking paper, 'C' was adjusted to match
    # experimental results for gallium melting.
    c_benchmark = 1.6 * 10**3

    # The question asks for the value from the "prior published implementation".
    requested_value = c_prior

    print("The Carman-Kozeny source term is a key component of the enthalpy-porosity model for simulating melting.")
    print("It takes the general form: Source = -C * [(1-f)^2 / (f^3 + b)] * u")
    print("where 'C' is the computational factor (mushy zone constant).")
    print("-" * 50)
    print(f"In the prior, simulation-only paper (Voller & Prakash, 1987), the value was C = {c_prior:.1e}.")
    print(f"In the later paper benchmarking against gallium melting (Brent, Voller & Reid, 1988), it was C = {c_benchmark:.1e}.")
    print("-" * 50)
    print("The question asks for the value from the prior published simulation-only work.")
    print("\nThe final answer is the value of 'C' from that original implementation.")
    print(f"\nC = {requested_value:.1e}")

# Execute the function to find and print the answer.
find_computational_factor()