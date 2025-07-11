def solve_cfgs():
    """
    This function analyzes three categories fibered in groupoids and prints their properties.
    """

    # Properties for X1: Hilb^11(A^3)
    # Type: Scheme (S), Separated (s), not uc, not irr, dim = 11*3=33
    x1_props = ['S', 's', 33]

    # Properties for X2: [ (A^4 \ V(xy-zw)) / C* ]
    # Type: Scheme (S), Separated (s), not uc, Irreducible (irr), dim = 4-1=3
    x2_props = ['S', 's', 'irr', 3]

    # Properties for X3: Pic(C_0) for a genus 7 curve C_0
    # Type: Algebraic Stack (A), Separated (s), not uc, not irr, dim = g = 7
    x3_props = ['A', 's', 7]

    # Helper function to format the properties list into the required string format
    def format_props(props):
        """Formats a list of properties into a bracketed, comma-separated string."""
        # This loop ensures each property, including the dimension, is included in the output string
        prop_strings = []
        for p in props:
            prop_strings.append(str(p))
        return f"[{','.join(prop_strings)}]"

    # Format each set of properties
    x1_str = format_props(x1_props)
    x2_str = format_props(x2_props)
    x3_str = format_props(x3_props)

    # Combine the formatted strings into the final result
    final_answer = f"{x1_str} {x2_str} {x3_str}"
    print(final_answer)

solve_cfgs()