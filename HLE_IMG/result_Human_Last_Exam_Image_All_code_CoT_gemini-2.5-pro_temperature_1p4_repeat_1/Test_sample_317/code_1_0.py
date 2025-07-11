import math

def solve_navier_stokes_puzzle():
    """
    This function analyzes the provided plots of the Navier-Stokes system
    to determine the 9-character solution string.
    The logic is explained in the comments.
    """

    # --- Part 1: Determine the Reynolds number constant 'k' ---
    # The equation for x3 is x3' = -5*x3 - c*x1*x2 + Re.
    # Taking the time-average (<.>) over a stable attractor gives: 0 = -5*<x3> - c*<x1*x2> + Re
    # Since Re = 50*k, this means <x3> = (50*k - c*<x1*x2>)/5 = 10*k - (c/5)*<x1*x2>.
    #
    # In Simulation 4, the vertical range of x3 (plots h, i) shrinks dramatically compared to Simulation 1.
    # The oscillating term in the x3' equation is -c*x1*x2. The x3 range in Sim 4 (~0.3) is about
    # 10 times smaller than in Sim 1 (~3.0). This strongly indicates the change in Sim 4 is a tenfold
    # decrease in parameter 'c', and the baseline value for c must be 10 (since initial parameters are -1 or 10).
    #
    # Now we can find 'k'. Let's use the average value of x3 from Sim 1, <x3> ≈ 16.5.
    # Test k=3 (Re=150): 16.5 = 10*3 - (10/5)*<x1*x2> => 16.5 = 30 - 2*<x1*x2> => <x1*x2> = 6.75.
    # This is impossible, as the maximum magnitude of x1*x2 in Sim 1 is |(-2)*(2)| = 4.
    # Test k=2 (Re=100): 16.5 = 10*2 - (10/5)*<x1*x2> => 16.5 = 20 - 2*<x1*x2> => <x1*x2> = 1.75.
    # This value is plausible given the ranges of x1 and x2 in Sim 1. Thus, k=2.
    k = 2

    # --- Part 2: Identify the plot axes ---
    # We established that x3 has a large positive average. The horizontal axis of plot 'f' is the
    # only one that is consistently large and positive across all simulations.
    # Therefore, the horizontal axis for x3 is 'f'.
    #
    # The other variables (x1, x2, x4) oscillate around zero. Their damping terms are -2*x1, -9*x2, -5*x4.
    # Weaker damping implies larger oscillations. Damping strength: x2 (9) > x4 (5) > x1 (2).
    # So, we expect amplitude order: |x1| > |x4| > |x2|.
    # In Sim 2, the ranges are distinct: horiz(h) is ~30, horiz(g) is ~12, horiz(i) is ~8.
    # This matches the expected order perfectly.
    # x1 -> h (largest amplitude), x4 -> g (middle), x2 -> i (smallest).
    axes = {'x1': 'h', 'x2': 'i', 'x3': 'f', 'x4': 'g'}
    axis_string = axes['x1'] + axes['x2'] + axes['x3'] + axes['x4']

    # --- Part 3: Identify the parameter changes ---
    # Sim 1 is the baseline with initial parameters.
    change1 = '0'
    # Sim 2: The ratio |x2|/|x1| is approx. |b|*|x3|/9.
    # In Sim 1, |x2|/|x1| ≈ 2/2 = 1. In Sim 2, |x2|/|x1| ≈ 4/15 ≈ 0.27.
    # Since the ratio decreased, |b| must have decreased. Change is 'b'.
    change2 = 'b'
    # Sim 3: x5 flips from negative to positive. This indicates the system jumped to a partner
    # attractor allowed by the (x4,x5) -> (-x4,-x5) symmetry. The change likely affects the (x4,x5)
    # subsystem (params d, e). The amplitude of x5 grows more than x4, pointing to a change in 'e'.
    # An increase is needed to drive the larger oscillations. Change is 'E'.
    change3 = 'E'
    # Sim 4: As established in Part 1, the x3 oscillation amplitude shrinks by a factor of 10.
    # This is caused by a tenfold decrease in parameter 'c'. Change is 'c'.
    change4 = 'c'

    change_string = change1 + change2 + change3 + change4

    # --- Part 4: Assemble the final string ---
    final_answer = str(k) + axis_string + change_string
    print(final_answer)

solve_navier_stokes_puzzle()