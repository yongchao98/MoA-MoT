import math

def explain_chemical_potential_limit():
    """
    Explains the fundamental limit on the chemical potential for bosons in Bose-Einstein condensation.
    """

    print("To find the limit on the chemical potential (μ) for bosons, we analyze the Bose-Einstein distribution.")
    print("The distribution gives the average number of particles, n(ε), in a state with energy ε:\n")
    print("  n(ε) = 1 / (exp((ε - μ) / (k_B * T)) - 1)\n")
    print("where T is temperature and k_B is Boltzmann's constant.")

    print("A core physical principle is that the number of particles n(ε) cannot be negative.")
    print("For n(ε) to be positive, the denominator of the fraction must be positive:\n")
    print("  exp((ε - μ) / (k_B * T)) - 1 > 0\n")

    print("This inequality is true only if the argument of the exponential function is positive:\n")
    print("  (ε - μ) / (k_B * T) > 0")
    print("  ==>  ε - μ > 0")
    print("  ==>  μ < ε\n")

    print("This condition, μ < ε, must hold for ALL possible energy states.")
    print("Therefore, it must hold for the lowest possible energy state, the ground state energy, ε_0.")
    print("This leads to the fundamental limit for the chemical potential of a Bose gas.")
    print("The final equation for this limit is:\n")
    
    # Define symbols for the final equation
    mu = "μ"
    less_than = "<"
    epsilon_0 = "ε_0"
    
    print(f"  {mu} {less_than} {epsilon_0}\n")
    
    print("This means the chemical potential must always be less than the ground state energy.")
    print("At absolute zero (T=0), a non-interacting Bose gas has all particles in the ground state,")
    print("and its chemical potential becomes exactly equal to the ground state energy, μ(T=0) = ε_0.")
    print("Thus, the ground state energy, which is the chemical potential of the gas at T=0, acts as the ultimate upper limit.")

explain_chemical_potential_limit()