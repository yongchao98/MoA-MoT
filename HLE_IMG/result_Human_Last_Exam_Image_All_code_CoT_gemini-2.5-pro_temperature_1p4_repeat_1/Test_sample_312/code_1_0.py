import sys

def solve_nanotube_puzzle():
    """
    This script determines the chiral index m for each of the nine plots
    and prints the final sequence as requested.
    The logic is based on analyzing the electronic properties of (4,m) SWNTs.

    - Plot 1: Metallic (m=1 or m=4), lower density of states -> m=1
    - Plot 2: Semiconducting (m=0,2,3), medium distance of cutting line to K-point -> m=2
    - Plot 3: Semiconducting (m=0,2,3), smallest band gap -> m=3
    - Plot 4: Metallic (cutting line on K-point), must be the m=1 case as the other metallic plot (m=4) is missing. -> m=1
    - Plot 5: Metallic (m=1 or m=4), higher density of states -> m=4
    - Plot 6: Semiconducting (m=0,2,3), smallest distance of cutting line to K-point -> m=3
    - Plot 7: Semiconducting (m=0,2,3), medium band gap -> m=2
    - Plot 8: Semiconducting (m=0,2,3), largest distance of cutting line to K-point -> m=0
    - Plot 9: Semiconducting (m=0,2,3), largest band gap -> m=0
    """
    
    # The determined sequence of m values for plots 1 through 9.
    m_values = [1, 2, 3, 1, 4, 3, 2, 0, 0]
    
    # Format the output as a sequence of nine integers in curly braces.
    # The puzzle requires outputting each number, so we join them.
    result_string = "{"
    for i, m in enumerate(m_values):
        result_string += str(m)
        if i < len(m_values) - 1:
            result_string += ", "
    result_string += "}"
    
    print(result_string)

solve_nanotube_puzzle()