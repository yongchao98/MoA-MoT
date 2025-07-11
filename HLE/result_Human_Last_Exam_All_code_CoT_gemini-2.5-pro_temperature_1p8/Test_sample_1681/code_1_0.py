import math

def solve_quantum_prisoners_dilemma():
    """
    Calculates the optimal equilibrium point for the quantum Prisoner's Dilemma.
    
    This function uses the established result from the Eisert-Wilkens-Lewenstein (EWL)
    quantum game theory framework.
    """
    
    # Step 1: Define the payoffs from the provided matrix.
    # The matrix is [[(R,R), (S,T)], [(T,S), (P,P)]]
    # R: Reward for mutual cooperation
    # S: Sucker's payoff if you cooperate and the other defects
    # T: Temptation to defect if the other cooperates
    # P: Punishment for mutual defection
    R = 5
    S = 0
    T = 7
    P = 1

    # Step 2: State the theoretical condition for the optimal equilibrium.
    # In the quantum game, (Cooperate, Cooperate) becomes a Nash Equilibrium if the
    # entanglement parameter γ satisfies the following condition:
    # tan(γ/2)^2 < (R - S) / (T - P)
    # This creates a game with two Nash Equilibria: (C,C) and (D,D).
    
    # Step 3: Calculate the critical value for the condition.
    threshold_numerator = R - S
    threshold_denominator = T - P
    
    # The condition is only meaningful if the denominator is not zero.
    if threshold_denominator == 0:
        print("The analysis cannot proceed as the Temptation and Punishment payoffs are equal.")
        return
        
    k_critical = threshold_numerator / threshold_denominator

    # Step 4: Explain the strategy based on the problem statement.
    # We are asked to "design the initial states for an optimal" outcome.
    # Choosing any level of entanglement 0 < γ < 2 * arctan(sqrt(k_critical))
    # ensures that (Cooperate, Cooperate) is a stable strategy.
    # Since the payoff for (C,C) is (R,R) = (5,5) and the payoff for (D,D) is
    # (P,P) = (1,1), rational players will coordinate on the superior (C,C) equilibrium.
    
    # Step 5: The resulting equilibrium point is the payoff for mutual cooperation.
    payoff_A = R
    payoff_B = R

    # Step 6: Print the detailed result and the final equation.
    print("In the quantum Prisoner's Dilemma, an optimal initial state can be designed.")
    print("The goal is to make the Pareto-optimal (Cooperate, Cooperate) outcome a Nash Equilibrium.")
    print("\nThis is achieved if the entanglement parameter γ satisfies the condition:")
    print(f"tan(γ/2)^2 < ({R} - {S}) / ({T} - {P})")
    print(f"tan(γ/2)^2 < {threshold_numerator} / {threshold_denominator}")
    print(f"tan(γ/2)^2 < {k_critical:.4f}")
    
    print("\nBy setting the initial entanglement within this range, rational players will choose to cooperate.")
    print("The equilibrium point is therefore the payoff vector for mutual cooperation.")
    print("\nFinal Equilibrium Payoff Equation:")
    print(f"Player A's Payoff = {payoff_A}")
    print(f"Player B's Payoff = {payoff_B}")

solve_quantum_prisoners_dilemma()