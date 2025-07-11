import math

def solve_factor_problem():
    """
    This function outlines the research to find the original computational factor
    used in the enthalpy-porosity method for simulating melting.
    """

    # The enthalpy-porosity method adds a source term to the momentum equation
    # to account for the mushy region, based on the Carman-Kozeny relations.
    # The source term S is often written as: S = -C * ((1 - f_l)^2 / (f_l^3 + b)) * u
    # where C is the computational factor, also called the mushy region constant.

    # Step 1: Identify the prior, simulation-only publication.
    # This is "A fixed grid numerical modelling methodology for convection-diffusion mushy region
    # phase-change problems" by V.R. Voller and C. Prakash (1987).
    prior_paper = "Voller and Prakash (1987)"

    # Step 2: Find the value of the factor 'C' in this paper.
    # In this paper, the authors state the value used for their general numerical model.
    original_computational_factor = 1.6 * (10**6)

    # Step 3: Identify the later, benchmarked publication.
    # This is "An enthalpy-porosity technique for modeling convection-diffusion phase change:
    # Application to the melting of a pure metal" by A.D. Brent, V.R. Voller, and K.J. Reid (1988).
    benchmark_paper = "Brent, Voller, and Reid (1988)"
    
    # Step 4: Find the modified value of 'C' in the benchmark paper.
    modified_computational_factor = 1.6 * (10**3)

    # The question asks for the value from the *prior* published work.
    final_answer = original_computational_factor
    
    print(f"The prior simulation-only paper, {prior_paper}, established the numerical method.")
    print(f"In this paper, the mushy region constant (computational factor) was set to {original_computational_factor:.1e}.")
    print("-" * 60)
    print(f"The later paper, {benchmark_paper}, benchmarked this method against gallium melting.")
    print(f"In this paper, the factor was modified to {modified_computational_factor:.1e} to match experiments.")
    print("-" * 60)
    print("The question asks for the original value from the prior implementation.")
    print("\nFinal Answer Equation:")

    base = 1.6
    exponent = int(math.log10(final_answer / base))
    
    print(f"Computational Factor = {base} * 10^{exponent}")

solve_factor_problem()