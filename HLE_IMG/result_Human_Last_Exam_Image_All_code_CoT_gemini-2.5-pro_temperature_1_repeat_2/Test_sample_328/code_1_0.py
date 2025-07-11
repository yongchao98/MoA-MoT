import collections

def solve_puzzle():
    """
    Solves the nonlinear wave equation parameter identification puzzle.

    The solution is determined by analyzing the visual characteristics of each plot
    and matching them to the physical effects of the parameters b, c, and d in the equation:
    d^2(Phi)/dt^2 = d^2(Phi)/dx^2 - Phi^3 + b*Phi^2 + c*Phi + d

    The reasoning for each plot's classification is as follows:

    - Parameter 'c' (stability):
      - c=-1 (code 'c'): Makes Phi=0 unstable, leading to chaotic, space-filling patterns.
        Plots #1, #2, #7, #11 show this characteristic.
      - c=1 (code 'C'): Makes Phi=0 stable, leading to regular, repeating oscillatory patterns.
        Plots #4, #6, #8, #9, #14, #16 show this characteristic.

    - Parameter 'd' (constant shift):
      - d=-1 (code 'd'): A constant negative force pushes the solution to negative (blue) values.
        Plots #5, #12, #13 clearly show a dominant blue background.
      - d=1 (code 'D'): A constant positive force pushes the solution to positive (yellow/red) values.
        Plot #15 has a distinct light yellowish background.

    - Parameter 'b' (asymmetry):
      - b=1 (code 'B'): Introduces an asymmetry favoring positive (red) values.
        Plot #3 shows the wave collision creating a prominent, persistent red structure.
      - b=-1 (code 'b'): Introduces an asymmetry favoring negative (blue) values.
        Plot #10 shows the wave collision creating a prominent blue V-shape.
        
    - The remaining codes (z, Z, 0) for unique zero-value parameters do not appear to
      correspond to any of the provided plots based on these primary visual cues.
    """

    # Mapping from plot number to identified code
    # The order of assignments follows the step-by-step reasoning.
    assignments = {
        # c=-1 (unstable, chaotic)
        1: 'c',
        2: 'c',
        7: 'c',
        11: 'c',
        # c=1 (stable, oscillatory)
        4: 'C',
        6: 'C',
        8: 'C',
        9: 'C',
        14: 'C',
        16: 'C',
        # d=-1 (strong negative/blue shift)
        5: 'd',
        12: 'd',
        13: 'd',
        # d=1 (positive/yellow shift)
        15: 'D',
        # b=1 (asymmetry favoring positive/red)
        3: 'B',
        # b=-1 (asymmetry favoring negative/blue)
        10: 'b'
    }

    # Sort the assignments by plot number to construct the final string
    sorted_assignments = collections.OrderedDict(sorted(assignments.items()))
    
    # Construct the final 16-character string
    result_string = "".join(sorted_assignments.values())
    
    print(result_string)

solve_puzzle()