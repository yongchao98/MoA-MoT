import numpy as np
import pandas as pd

# Set a seed for reproducibility
np.random.seed(42)

# Number of data points
n = 10000

# --- 1. Define the structural equations and generate data ---
# E is the common cause
E = np.random.normal(loc=0, scale=2, size=n)

# A is caused by E
# E -> A
noise_A = np.random.normal(loc=0, scale=1, size=n)
A = 0.8 * E + noise_A

# D is also caused by E
# E -> D
noise_D = np.random.normal(loc=0, scale=1, size=n)
D = 0.7 * E + noise_D

# B is caused by A
# A -> B
noise_B = np.random.normal(loc=0, scale=1, size=n)
B = 0.6 * A + noise_B

# C is caused by B and D (C is a collider)
# B -> C <- D
noise_C = np.random.normal(loc=0, scale=1, size=n)
C = 0.5 * B + 0.4 * D + noise_C


# --- 2. Create a DataFrame and calculate correlations ---
df = pd.DataFrame({'A': A, 'B': B, 'C': C, 'D': D, 'E': E})
correlation_matrix = df.corr()

# Get the correlation between A and D
corr_AD = correlation_matrix.loc['A', 'D']


# --- 3. Print the results and the conclusion ---
print("Simulated Structural Equation Model: E->A->B->C<-D<-E")
print("\nIn our simulation, we defined the causal links as follows:")
print("1. A is caused by E (A = 0.8*E + noise)")
print("2. D is caused by E (D = 0.7*E + noise)")
print("3. There is NO direct or indirect causal path between A and D.")

print(f"\nThe calculated correlation between A and D is: {corr_AD:.4f}")

print("\nAs you can see, the correlation is significantly different from zero.")
print("This correlation is not due to a causal link between A and D.")
print("It is a spurious correlation caused by the common cause 'E'.")

print("\nTherefore, does correlation imply causation in this system?")
print("The answer is:")

<<<No>>>