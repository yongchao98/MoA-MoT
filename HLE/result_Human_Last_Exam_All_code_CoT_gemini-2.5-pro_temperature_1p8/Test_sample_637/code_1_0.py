import numpy as np

def design_deadbeat_observer():
    """
    Designs a deadbeat observer for the given system.

    For a discrete-time system x(k + 1) = Ax(k) + Bu(k), y(k) = Cx(k),
    the observer dynamics are given by x_hat(k+1) = A*x_hat(k) + B*u(k) + L*(y(k) - C*x_hat(k)).
    The error dynamics are e(k+1) = (A - LC)e(k).

    A deadbeat observer requires the eigenvalues of (A - LC) to be all zero.
    This makes the characteristic polynomial of (A - LC) equal to z^n = 0.

    We will solve for the observer gain L. We assume a simplified L that only uses
    the first output channel, i.e., the second column of L is zero.
    L = [[l1, 0], [l2, 0], [l3, 0], [l4, 0]]

    The characteristic polynomial of (A - LC) for this form of L was derived as:
    p(z) = z^4 + (l1+2)z^3 + (l1+l4+3)z^2 + (l1+l3+3)z + (3+2*l1-l2)

    To make p(z) = z^4, we set the coefficients of the lower-order terms to zero:
    1. Coeff of z^3: l1 + 2 = 0
    2. Coeff of z^2: l1 + l4 + 3 = 0
    3. Coeff of z^1: l1 + l3 + 3 = 0
    4. Coeff of z^0: 3 + 2*l1 - l2 = 0

    We solve this system of linear equations for l1, l2, l3, l4.
    """

    # Create a vector for the elements of the first column of L
    l = np.zeros(4)

    # From equation 1:
    l[0] = -2  # l1

    # From equation 2:
    l[3] = -l[0] - 3  # l4

    # From equation 3:
    l[2] = -l[0] - 3  # l3

    # From equation 4:
    l[1] = 3 + 2 * l[0]  # l2

    # Construct the observer gain matrix L
    L = np.zeros((4, 2))
    L[:, 0] = l

    print("The observer gain matrix L is:")
    print(L)

design_deadbeat_observer()