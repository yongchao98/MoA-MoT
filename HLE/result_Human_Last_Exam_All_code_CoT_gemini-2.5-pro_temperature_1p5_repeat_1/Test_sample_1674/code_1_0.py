import numpy as np

def run_demonstration():
    """
    Demonstrates why adding a polarizing birefringent medium can make an optical
    system non-invertible, causing the time-reversal theory to fail.
    """
    print("This script demonstrates the failure of the inversion theory when a polarizing medium is added.")
    print("-" * 70)

    # A birefringent medium can act as a polarizer.
    # We model a perfect horizontal polarizer with a Jones matrix.
    # This matrix projects any input polarization onto the horizontal axis.
    T_biref = np.array([[1, 0],
                        [0, 0]])

    # Let's define a simple, invertible matrix for the initial "random medium".
    # This could be any transformation that doesn't lose information.
    T_rand = np.array([[0.8, 0.2],
                       [-0.2, 0.8]])

    # The total transmission matrix of the medium is the combination of the two.
    # The order doesn't matter for demonstrating non-invertibility.
    T_total = T_biref @ T_rand

    print("The Jones matrix for the added polarizing birefringent layer is:")
    print(T_total)
    print("\nLet's check the numbers in this matrix equation:")
    print(f"T_total = [[1, 0], [0, 0]] @ [[{T_rand[0,0]}, {T_rand[0,1]}], [{T_rand[1,0]}, {T_rand[1,1]}]]")
    print(f"T_total = [[{T_total[0,0]}, {T_total[0,1]}], [{T_total[1,0]}, {T_total[1,1]}]]")
    print("-" * 70)


    # The theory of inverting the system depends on the invertibility of T_total.
    # A matrix is invertible if and only if its determinant is non-zero.
    determinant = np.linalg.det(T_total)

    print(f"The determinant of the total transmission matrix is: {determinant}")
    print("\nFinal equation for the determinant:")
    print(f"det(T_total) = ({T_total[0,0]} * {T_total[1,1]}) - ({T_total[0,1]} * {T_total[1,0]})")
    print(f"det(T_total) = {determinant}")
    print("-" * 70)

    if determinant == 0:
        print("Since the determinant is zero, the matrix is singular (non-invertible).")
        print("This means information is lost permanently.")
    else:
        print("The matrix is invertible.")

    # Now, let's try to find the inverse matrix needed for the theory to hold.
    print("\nAttempting to calculate the inverse of the total transmission matrix...")
    try:
        T_inv = np.linalg.inv(T_total)
        print("Inverse matrix was found (this should not happen in this demo):")
        print(T_inv)
    except np.linalg.LinAlgError as e:
        print(f"Failed to compute the inverse. Reason: {e}")
        print("\nConclusion: Because the inverse matrix does not exist, the optical system cannot be inverted.")
        print("Therefore, the theory that you can recover the necessary input from the output does not hold.")

if __name__ == '__main__':
    run_demonstration()
