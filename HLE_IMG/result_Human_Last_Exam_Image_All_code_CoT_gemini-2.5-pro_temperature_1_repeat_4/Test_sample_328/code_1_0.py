def solve_puzzle():
    """
    Solves the IPhO 2011 nonlinear wave equation problem.

    The solution is determined by analyzing the effects of parameters b, c, and d
    on the nonlinear wave equation and matching them to the 16 plots. The logic is as follows:

    1.  Plot 1: (b=1, c=0, d=0). Unique b=1 -> 'B'
    2.  Plot 2: (b=-1, c=0, d=0). Unique b=-1 -> 'b'
    3.  Plot 3: (d=1, b=0, c=0). Unique d=1 -> 'D'
    4.  Plot 4: (c=-1, b=0, d=0). Unique c=-1 -> 'c'. This is a symmetric, stable case.
    5.  Plot 5: (b=0, c=-1, d=-1). Unique b=0 -> 'z'. Background decay to blue.
    6.  Plot 6: (b=-1, c=1, d=1). Unique b=-1 -> 'b'.
    7.  Plot 7: (b=0, c=1, d=1). Unique b=0 -> 'z'.
    8.  Plot 8: (c=0, b=1, d=1). Unique c=0 -> 'Z'.
    9.  Plot 9: (c=0, b=-1, d=-1). Unique c=0 -> 'Z'.
    10. Plot 10: (c=1, b=0, d=0). Unique c=1 -> 'C'. Classic kink-formation.
    11. Plot 11: (b=1, c=-1, d=-1). Unique b=1 -> 'B'.
    12. Plot 12: (d=-1, b=0, c=0). Unique d=-1 -> 'd'. Clean background decay to blue.
    13. Plot 13: (c=1, b=-1, d=-1). Unique c=1 -> 'C'.
    14. Plot 14: (d=0, b=1, c=1). Unique d=0 -> '0'.
    15. Plot 15: (d=0, b=-1, c=-1). Unique d=0 -> '0'.
    16. Plot 16: (c=-1, b=1, d=1). Unique c=-1 -> 'c'.

    This mapping is based on a detailed analysis of the equation's symmetries,
    equilibria, and stability, cross-referenced with the visual evidence in each plot.
    """
    
    # The final answer string is constructed based on the identification for each plot.
    answer_string = "BbDczbzZZCBdC00c"
    print(answer_string)

solve_puzzle()