def solve_reaction_diffusion_plots():
    """
    Analyzes the six plots and determines the parameter change for each.

    The logic is as follows:
    1.  Identify Plot #2 as the baseline (0) due to its intermediate characteristics.
    2.  Identify changes in the Diffusion coefficient (D):
        - Plot #5 shows a very flat concentration profile for A, indicating diffusion is much faster than reaction. This is D (doubled D).
        - Plot #6 shows a very steep profile for A and low product B, indicating A cannot diffuse to the center. This is d (halved D).
    3.  Identify changes in reaction parameters (k, n):
        - Plot #3 shows the most extreme reaction, with A nearly depleted and a sharp "knee" in the profile, characteristic of a high reaction order. This is N (doubled n).
        - Plot #1 shows a faster reaction than the baseline (lower A, higher B). This is K (doubled k).
        - Plot #4 shows a slower reaction than the baseline (higher A, lower B). This is k (halved k).
    4.  The parameter change 'n' (halved reaction order) is not represented in the plots.
    5.  Assemble the string by ordering the codes for plots 1 through 6.
    """
    # Plot 1: Faster reaction than baseline -> K
    # Plot 2: Baseline -> 0
    # Plot 3: Very fast, sharp-knee profile -> N
    # Plot 4: Slower reaction than baseline -> k
    # Plot 5: Very flat profile (fast diffusion) -> D
    # Plot 6: Very steep profile (slow diffusion) -> d
    
    answer = "K0NkDd"
    print(answer)

solve_reaction_diffusion_plots()