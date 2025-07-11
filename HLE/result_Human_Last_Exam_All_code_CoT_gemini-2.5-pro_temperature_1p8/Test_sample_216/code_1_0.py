import math

def solve_imitation_learning_bound():
    """
    This function calculates and prints the symbolic upper bound for the
    performance difference in an imitation learning setting based on the provided inputs.
    """

    # Symbolic representations of the variables
    H = 'H'
    A_size = '|A|'
    lmbda = 'λ'

    # The logic is based on the performance difference bound in imitation learning:
    # J(pi*) - J(pi_hat) <= H * R_max * D_TV(d_pi*, d_pi_hat).
    # Assuming R_max = 1 (rewards normalized to [0,1]).
    # The "population TV risk" T(pi_hat, pi*) is interpreted as D_TV(d_pi*, d_pi_hat).
    # The problem gives T(pi_hat, pi*) <= |A|*(1 - exp(-λ)).
    # Combining these gives the final bound.

    print("The tightest upper bound of J(pi^*) - J(pi_hat) is expressed by the formula:")

    # Printing each symbolic part of the final equation
    print(f"({H}) * ({A_size}) * (1 - exp(-{lmbda}))")

if __name__ == '__main__':
    solve_imitation_learning_bound()