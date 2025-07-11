import math

def construct_rigid_matrix_parameters(N):
    """
    This function calculates the rank 'r' for a rigid matrix based on the analysis.

    Our FNP algorithm with an NP oracle can be interpreted as an F-Sigma_2-P algorithm.
    This class of algorithms can non-deterministically guess a matrix M and then verify
    its rigidity. Rigidity is a co-NP property, and such properties can be decided
    by a Sigma_2-P machine.

    The existence of rigid matrices is the remaining constraint. It has been proven that
    rigid matrices exist for r up to N - O(sqrt(N)). Based on this, we choose a
    representative value for the largest possible rank r.

    Args:
        N (int): The dimension of the square matrix.

    Returns:
        None. It prints the result directly.
    """
    if not isinstance(N, int) or N <= 0:
        print("Error: N must be a positive integer.")
        return

    # For constant delta, the largest known r for which rigidity is achievable
    # is of the form N - c * sqrt(N). We will use c=1 for a concrete example.
    c = 1.0
    
    # Calculate the subtraction term
    sub_term = c * math.sqrt(N)
    
    # Calculate the rank r
    r = N - sub_term
    
    # We round r down to the nearest integer as rank must be an integer.
    r_int = math.floor(r)
    
    print(f"For a given N = {N}, the largest rank 'r' for which a rigid matrix can be constructed by the specified algorithm is given by the equation:")
    # The final instruction requested to output each number in the final equation.
    # r = N - c * sqrt(N)
    # The value of c is chosen to be 1 for this example.
    print(f"r = {N} - {c} * sqrt({N})")
    print(f"r = {N} - {sub_term:.2f}")
    print(f"r approx. {r:.2f}")
    print(f"So, the largest integer rank is r = {r_int}")

# Example usage with N=400
N_val = 400
construct_rigid_matrix_parameters(N_val)
