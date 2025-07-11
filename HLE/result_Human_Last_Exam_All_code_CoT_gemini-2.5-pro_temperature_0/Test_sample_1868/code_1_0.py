import sys

# Given parameters
c = 0.95  # Consistency level of identifier resolution
b = 3     # Branching factor of semantic version control

# --- Step 1: Model the FAIR metrics based on the plan ---

# Findability (f) is modeled as the consistency level c.
f = c

# Accessibility (a) is assumed to be ideal (1.0) under best practices.
a = 1.0

# Interoperability (i) is modeled as inversely proportional to the branching factor b.
# A higher branching factor complicates integration, reducing interoperability.
i = 1 / b

# Reusability (r) is assumed to be ideal (1.0) under best practices.
r = 1.0

# --- Step 2: Model the overall FAIR score R ---
# R is the arithmetic mean of the four metrics.
R = (f + a + i + r) / 4

# --- Step 3: Print the calculation and result ---
print("Calculating the FAIR Compliance Score (R)")
print("-----------------------------------------")
print(f"Model: R = (Findability + Accessibility + Interoperability + Reusability) / 4")
print(f"Findability (f) = c = {f:.2f}")
print(f"Accessibility (a) = Best Practice = {a:.2f}")
print(f"Interoperability (i) = 1/b = 1/{b} = {i:.4f}")
print(f"Reusability (r) = Best Practice = {r:.2f}")
print("\nFinal Equation:")
# The final print statement shows the full equation with the calculated numbers
print(f"R = ({f:.2f} + {a:.2f} + {i:.4f} + {r:.2f}) / 4 = {R:.4f}")

# The problem asks for the final answer in a specific format.
# To avoid printing this line in the user's terminal, we can write it to stderr
# or just accept that it will be part of the script's output.
# For this task, we will just print it.
sys.stdout.write(f"\n<<<{R:.4f}>>>\n")