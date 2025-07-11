import numpy as np

def calculate_moments(E, m2):
    """
    Calculates even moments <x^{2k}> from k=0 to 7 using the recursion relation.
    
    Args:
        E (float): The test energy value.
        m2 (float): The test value for <x^2>.
        
    Returns:
        numpy.ndarray: An array containing the moments m_0, m_2, ..., m_14.
    """
    # We need moments up to m_14, so an array of size 15 for m_0 to m_14 is sufficient.
    # Odd moments are zero and not stored.
    m = np.zeros(16) 
    m[0] = 1.0
    m[2] = m2

    # The general recursion relation is:
    # (4t+8) m_{t+3} = 4tE m_{t-1} - (4t+4) m_{t+1} + t(t-1)(t-2) m_{t-3}
    # We use it for odd t to get relations between even moments.

    # t=1: (12) m_4 = 4*E*m_0 - 8*m_2
    m[4] = (4 * E * m[0] - 8 * m[2]) / 12

    # t=3: (20) m_6 = 12*E*m_2 - 16*m_4 + 6*m_0
    m[6] = (12 * E * m[2] - 16 * m[4] + 6 * m[0]) / 20

    # t=5: (28) m_8 = 20*E*m_4 - 24*m_6 + 60*m_2
    m[8] = (20 * E * m[4] - 24 * m[6] + 60 * m[2]) / 28

    # t=7: (36) m_{10} = 28*E*m_6 - 32*m_8 + 210*m_4
    m[10] = (28 * E * m[6] - 32 * m[8] + 210 * m[4]) / 36

    # t=9: (44) m_{12} = 36*E*m_8 - 40*m_{10} + 504*m_6
    m[12] = (36 * E * m[8] - 40 * m[10] + 504 * m[6]) / 44

    # t=11: (52) m_{14} = 44*E*m_{10} - 48*m_{12} + 990*m_8
    m[14] = (44 * E * m[10] - 48 * m[12] + 990 * m[8]) / 52
    
    return m

def check_psd(moments):
    """
    Checks if the moment matrices M_even and M_odd are positive semidefinite.
    
    Args:
        moments (numpy.ndarray): An array of even moments.
        
    Returns:
        bool: True if both matrices are positive semidefinite, False otherwise.
    """
    # For K=7, the matrices are 4x4.
    # M_even entries are <x^{2(i+j)}> = m_{2(i+j)}
    # M_odd entries are <x^{2(i+j+1)}> = m_{2(i+j+1)}
    
    m_even = np.zeros((4, 4))
    m_odd = np.zeros((4, 4))
    
    for i in range(4):
        for j in range(4):
            m_even[i, j] = moments[2 * (i + j)]
            m_odd[i, j] = moments[2 * (i + j + 1)]

    # A matrix is positive semidefinite if all its eigenvalues are non-negative.
    # We use a small tolerance for floating point inaccuracies.
    try:
        eigvals_even = np.linalg.eigvalsh(m_even)
        eigvals_odd = np.linalg.eigvalsh(m_odd)
    except np.linalg.LinAlgError:
        return False

    return np.all(eigvals_even >= -1e-9) and np.all(eigvals_odd >= -1e-9)

def find_ground_state():
    """
    Performs a grid search to find the minimal E and corresponding <x^2>.
    """
    # Define the search grid with a step of 0.001 for 3-digit precision.
    e_range = np.arange(1.3, 1.5, 0.001)
    m2_range = np.arange(0.3, 0.5, 0.001)

    # Iterate from low E to high E. The first valid point found will be the minimum.
    for E in e_range:
        for m2 in m2_range:
            moments = calculate_moments(E, m2)
            if check_psd(moments):
                # This is the first point in the allowed region found by the search.
                # It corresponds to the ground state energy and expectation value.
                print(f"Minimal E = {E:.3f}")
                print(f"Minimal <x^2> = {m2:.3f}")
                print(f"<<<E={E:.3f}, <x^2>={m2:.3f}>>>")
                return

    print("Solution not found in the specified search range.")

if __name__ == '__main__':
    find_ground_state()