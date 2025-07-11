import math

def nCr(n, r):
    """Calculates the binomial coefficient 'n choose r'."""
    f = math.factorial
    return f(n) // (f(r) * f(n - r))

# --- Introduction ---
print("Calculating symmetry factors for second-order vacuum bubbles in phi^4 theory.")
print("The theory has a L_int = -lambda/4! * phi^4 interaction.")
print("Second-order diagrams contain V=2 vertices, each with v=4 legs.")
print("The symmetry factor S is calculated as S = (V! * (v!)^V) / N_wick, where N_wick is the number of Wick contractions for the diagram.")
print("-" * 40)

# --- Diagram 1: Setting-Sun Diagram ---
print("Diagram 1: The 'Setting-Sun' Diagram")
print("This diagram consists of two vertices connected directly by four propagators.")
print("Let's calculate its Wick contraction count (N_A):")
print("This corresponds to the number of ways to pair the 4 legs of vertex 1 with the 4 legs of vertex 2.")

# N_A is 4!
N_A = math.factorial(4)
print(f"The number of ways is 4! = {N_A}.")

# Calculate S_A
V = 2
v = 4
S_A = (math.factorial(V) * math.factorial(v)**V) / N_A
print(f"The symmetry factor S_A = (2! * (4!)^2) / {N_A} = (2 * {math.factorial(v)}^2) / {N_A} = {int(S_A)}.")
print("-" * 40)

# --- Diagram 2: Figure-Eight Diagram ---
print("Diagram 2: The 'Figure-Eight' Diagram")
print("This diagram consists of two vertices, each with one loop, connected by two propagators.")
print("Let's calculate its Wick contraction count (N_B):")

# Step 1: Choose legs for loop at vertex 1
ways_v1_loop = nCr(4, 2)
print(f"1. Choose 2 of 4 legs at vertex 1 to form a loop: C(4, 2) = {ways_v1_loop} ways.")

# Step 2: Choose legs for loop at vertex 2
ways_v2_loop = nCr(4, 2)
print(f"2. Choose 2 of 4 legs at vertex 2 to form a loop: C(4, 2) = {ways_v2_loop} ways.")

# Step 3: Connect the remaining legs
ways_connect = 2
print("3. Connect the 2 remaining legs from vertex 1 to the 2 at vertex 2: 2 ways.")

N_B = ways_v1_loop * ways_v2_loop * ways_connect
print(f"The total number of Wick contractions N_B = {ways_v1_loop} * {ways_v2_loop} * {ways_connect} = {N_B}.")

# Calculate S_B
S_B = (math.factorial(V) * math.factorial(v)**V) / N_B
print(f"The symmetry factor S_B = (2! * (4!)^2) / {N_B} = (2 * {math.factorial(v)}^2) / {N_B} = {int(S_B)}.")
print("-" * 40)

# --- Summation ---
print("Finally, we sum the symmetry factors of the two diagrams:")
total_S = S_A + S_B
print(f"Total Symmetry Factor = S_A + S_B = {int(S_A)} + {int(S_B)} = {int(total_S)}")
