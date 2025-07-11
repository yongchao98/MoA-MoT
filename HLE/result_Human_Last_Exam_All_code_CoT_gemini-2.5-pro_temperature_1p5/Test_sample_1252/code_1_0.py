def solve_moduli_question(g, num_markings):
    """
    Calculates and prints the properties of the moduli space M_trop(g, n).

    Args:
        g (int): The genus of the curve.
        num_markings (int): The number of marked points, |A|.
    """
    # Check if the moduli space is non-empty. The properties are typically
    # discussed for non-empty spaces, where 2g - 2 + n > 0.
    if 2 * g - 2 + num_markings <= 0:
        print(f"For g={g} and |A|={num_markings}, the moduli space M_trop(g,A) is empty or trivial.")
        return

    # Part (a): Expression for the number of vertices.
    # This is the number of vertices for a generic (trivalent) graph type.
    num_vertices = 2 * g - 2 + num_markings

    # Part (b): The answer is static.
    is_simplicial_fan = "yes"

    # Part (c): The answer depends on the dimension formula.
    is_tropical_variety = "yes"
    complex_dimension = 3 * g - 3 + num_markings

    # Format the final output string with the computed numbers.
    # The prompt requests each number in the final equation to be outputted.
    final_answer = (
        f"(a) {num_vertices}; "
        f"(b) {is_simplicial_fan}; "
        f"(c) {is_tropical_variety}, {complex_dimension}."
    )
    print(final_answer)

# --- Example Usage ---
# Let's solve for a genus 2 curve with 3 marked points.
print("Example for g=2, |A|=3:")
solve_moduli_question(g=2, num_markings=3)

# Let's solve for a genus 0 curve with 5 marked points.
print("\nExample for g=0, |A|=5:")
solve_moduli_question(g=0, num_markings=5)

# Let's solve for a genus 1 curve with 1 marked point.
print("\nExample for g=1, |A|=1:")
solve_moduli_question(g=1, num_markings=1)