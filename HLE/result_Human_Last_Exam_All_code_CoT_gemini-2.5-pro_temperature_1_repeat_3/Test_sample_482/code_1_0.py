import sys

def solve_critical_correlation():
    """
    This function derives and prints the formula for the critical amount of correlation
    between input populations v and s that balances potentiation and depression.

    The derivation is based on the following steps:
    1. The total expected change in synaptic weights for an output neuron must be zero.
       E[d(W_total)/dt] = 0.
    2. The dynamics are simplified by assuming symmetric input populations, i.e.,
       N_v = N_s = N and sigma_v^2 = sigma_s^2 = sigma_sq.
    3. This leads to the condition: sigma_sq + N*C_vs + 2*N*(mu^2 - theta*mu) = 0.
    4. Solving for C_vs gives the critical correlation.
    """

    # Define the symbols used in the equation
    # C_vs: The critical covariance between input populations v and s.
    # mu: The average firing rate of input neurons.
    # theta: The heterosynaptic depression constant.
    # sigma_sq: The variance of the firing rate of input neurons.
    # N: The number of neurons in each input population.

    print("The equation for the critical amount of correlation (C_vs) is:")
    
    # Use the print function to output each part of the final equation.
    # We use `end=''` to print everything on a single line.
    print("C_vs", end="")
    sys.stdout.flush()
    print(" = ", end="")
    sys.stdout.flush()
    print("2", end="")
    sys.stdout.flush()
    print(" * ", end="")
    sys.stdout.flush()
    print("mu", end="")
    sys.stdout.flush()
    print(" * ", end="")
    sys.stdout.flush()
    print("(theta - mu)", end="")
    sys.stdout.flush()
    print(" - ", end="")
    sys.stdout.flush()
    print("sigma_sq", end="")
    sys.stdout.flush()
    print(" / ", end="")
    sys.stdout.flush()
    print("N")
    sys.stdout.flush()

solve_critical_correlation()