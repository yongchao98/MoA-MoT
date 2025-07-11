def compute_edeg_bounds(variety_name, properties):
    """
    Computes the pair (m(X), M(X)) for a given variety X.
    This function encodes known results from algebraic geometry regarding the
    effective degree of 0-cycles on specific types of varieties.

    Args:
        variety_name (str): The name of the variety (e.g., 'X_1').
        properties (dict): A dictionary containing key properties of the variety.

    Returns:
        tuple: A tuple (m, M) for the variety.
    """
    variety_type = properties.get('type')

    if variety_type == 'curve':
        g = properties.get('genus')
        if g is None:
            raise ValueError("Genus must be specified for a curve.")

        # A genus 2 curve is always hyperelliptic.
        is_hyperelliptic = (g == 2)
        is_general = properties.get('is_general', False)

        if g >= 3 and is_general:
            # For a general curve C of genus g >= 3, edeg(C, p) = g+1 for all points p.
            m = g + 1
            M = g + 1
        elif is_hyperelliptic:
            # For a hyperelliptic curve of genus g, edeg is 2 at Weierstrass points
            # and g+1 at all other points.
            m = 2
            M = g + 1
        else:
            # This case is more complex and not needed for this problem.
            return ('unknown', 'unknown')
        
        return (m, M)

    elif variety_type == 'Enriques surface':
        # For an Enriques surface S, 2[p] is rationally equivalent to 2[q] for all p, q.
        # This implies edeg(S, p) = 2 for all points p.
        return (2, 2)

    elif variety_type == 'Grassmannian':
        # Grassmannians are rational varieties. For any rational variety X,
        # CH_0(X)_0 is trivial, which implies edeg(X, p) = 1 for all points p.
        return (1, 1)

    else:
        return ('unknown', 'unknown')

def solve_problem():
    """
    Solves the problem for the four given varieties and prints the result.
    """
    varieties = {
        'X_1': {'type': 'curve', 'genus': 2},
        'X_2': {'type': 'curve', 'genus': 7, 'is_general': True},
        'X_3': {'type': 'Enriques surface'},
        'X_4': {'type': 'Grassmannian'}
    }

    results = []
    for name in ['X_1', 'X_2', 'X_3', 'X_4']:
        properties = varieties[name]
        m_val, M_val = compute_edeg_bounds(name, properties)
        results.append(f"({m_val}, {M_val})")

    # The requested output format is "(m(X_1), M(X_1)), (m(X_2), M(X_2)), ..."
    final_output = ", ".join(results)
    print(final_output)

if __name__ == "__main__":
    solve_problem()