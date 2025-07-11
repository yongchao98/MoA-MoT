import math

def get_nth_moment_g(n):
    """
    Calculates the n-th moment of the dimensionless conductance 'g' for a class D
    disordered wire at the critical point.
    The universal formula is <g^n> = 1 / (2n + 1).
    """
    if n <= 0:
        raise ValueError("Moment number 'n' must be a positive integer.")
    return 1.0 / (2 * n + 1)

def main():
    """
    Calculates the ratio of the fourth moment of conductance to its average value.
    """
    # Calculate the average value (1st moment) of conductance
    n1 = 1
    avg_g = get_nth_moment_g(n1)

    # Calculate the 4th moment of conductance
    n4 = 4
    fourth_moment_g = get_nth_moment_g(n4)

    # Calculate the ratio
    ratio = fourth_moment_g / avg_g

    print("This script calculates the ratio of the fourth statistical moment of dimensionless conductance to its average value for a critical disordered Majorana wire.")
    print("-" * 50)
    print(f"The average conductance (1st moment, <g>) is: 1 / (2*{n1} + 1) = {avg_g:.4f}")
    print(f"The 4th moment of conductance (<g^4>) is: 1 / (2*{n4} + 1) = {fourth_moment_g:.4f}")
    print("-" * 50)
    print(f"The ratio <g^4> / <g> is: {fourth_moment_g:.4f} / {avg_g:.4f} = {ratio:.4f}")

if __name__ == "__main__":
    main()