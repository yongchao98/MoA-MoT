def find_symmetric_psne(n, incentive_function):
    """
    Finds all pure strategy Nash equilibria in an n-player symmetric game.

    Args:
        n (int): The number of players.
        incentive_function (list): A list of n numbers representing f(j) = u(1,j) - u(0,j)
                                   for j = 0, 1, ..., n-1.
    """
    print(f"Analyzing a {n}-player game...")
    print(f"Incentive function f(j) for j=0 to {n-1}: {incentive_function}\n")
    
    f = incentive_function
    psne_list = []

    # Iterate through all possible symmetric strategy profiles
    # k is the number of players choosing Action 1
    for k in range(n + 1):
        is_psne = False
        # Case 1: k=0 (all players choose Action 0)
        # Condition: No player wants to switch to Action 1.
        # This means the incentive f(0) must be non-positive.
        if k == 0:
            if f[0] <= 0:
                is_psne = True
                print(f"Checking k={k}: All players choose Action 0.")
                print(f"  - Condition: f(0) <= 0. We have f(0) = {f[0]}. Condition MET.")

        # Case 2: k=n (all players choose Action 1)
        # Condition: No player wants to switch to Action 0.
        # This means the incentive f(n-1) must be non-negative.
        elif k == n:
            if f[n-1] >= 0:
                is_psne = True
                print(f"Checking k={k}: All players choose Action 1.")
                print(f"  - Condition: f(n-1) >= 0. We have f({n-1}) = {f[n-1]}. Condition MET.")

        # Case 3: 0 < k < n (a mix of actions)
        # Two conditions must hold:
        # 1. Players choosing Action 1 don't want to switch: f(k-1) >= 0
        # 2. Players choosing Action 0 don't want to switch: f(k) <= 0
        else:
            condition1 = f[k-1] >= 0
            condition2 = f[k] <= 0
            if condition1 and condition2:
                is_psne = True
                print(f"Checking k={k}: {k} players choose Action 1, {n-k} choose Action 0.")
                print(f"  - Condition 1: f(k-1) >= 0. We have f({k-1}) = {f[k-1]}. Condition MET.")
                print(f"  - Condition 2: f(k) <= 0. We have f({k}) = {f[k]}. Condition MET.")

        if is_psne:
            psne_list.append(k)
            print(f"  -> RESULT: k={k} is a Pure Strategy Nash Equilibrium.\n")

    print("---" * 10)
    if not psne_list:
        print("No Pure Strategy Nash Equilibria found.")
    else:
        print(f"Found {len(psne_list)} PSNE profile(s) where k players choose Action 1:")
        print(f"k = {psne_list}")
    print("---" * 10)


if __name__ == "__main__":
    # Let's define a 5-player game (n=5)
    num_players = 5
    
    # We need to define the incentive function f(j) for j = 0, 1, 2, 3, 4.
    # Let's create a function that "crosses" from positive to negative only once
    # to create a single PSNE.
    # f = [10, 5, -2, -8, -15] is a good example.
    # f(0)=10, f(1)=5, f(2)=-2, f(3)=-8, f(4)=-15
    incentives = [10, 5, -2, -8, -15]
    
    find_symmetric_psne(num_players, incentives)
