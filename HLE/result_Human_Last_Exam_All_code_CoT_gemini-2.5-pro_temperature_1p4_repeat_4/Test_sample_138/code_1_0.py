# This code must be run in a SageMath environment,
# as it relies on its advanced number theory library.

from sage.lattices.root_lattice import RootLattice
from sage.quadratic_forms.integer_quadratic_form import IntegerQF

def solve_lattice_problem():
    """
    Finds the number of positive definite even lattices of dimension 17
    and determinant 2.
    """
    dimension = 17
    determinant = 2

    # Step 1: Construct a representative lattice in the desired genus.
    # We use the direct sum of root lattices: E8 + E8 + A1.
    # E8 is dimension 8, det 1, even.
    # A1 is dimension 1, det 2, even.
    # The sum E8 + E8 + A1 has:
    # dimension = 8 + 8 + 1 = 17
    # determinant = 1 * 1 * 2 = 2
    # It is also positive definite and even.
    
    # Get the Gram matrices from Sage's RootLattice factory
    E8_matrix = RootLattice("E8").gram_matrix()
    A1_matrix = RootLattice("A1").gram_matrix()

    # The Gram matrix of the direct sum is the block diagonal matrix
    from sage.matrix.constructor import block_diagonal_matrix
    gram_matrix = block_diagonal_matrix(E8_matrix, E8_matrix, A1_matrix)

    # Create the integer quadratic form object from the Gram matrix
    qf = IntegerQF(gram_matrix)

    # Step 2: Get the genus of this quadratic form.
    # This object represents all lattices locally isometric to ours.
    genus = qf.genus()

    # Step 3: Enumerate all distinct isometry classes in the genus.
    # For positive definite forms, Sage can compute all of them.
    isometry_classes = genus.representatives()

    # The number of classes is the answer.
    num_classes = len(isometry_classes)

    # Final output message
    print(f"The problem is to find the number of positive definite even lattices of dimension {dimension} and determinant {determinant}.")
    print(f"The final count is: {num_classes}")

solve_lattice_problem()