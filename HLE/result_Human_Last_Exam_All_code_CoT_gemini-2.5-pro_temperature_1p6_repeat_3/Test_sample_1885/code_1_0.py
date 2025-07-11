import random

# Finite stand-ins for transfinite cardinals
# W1 represents the domain/codomain size (like omega_1)
W1 = 10
# W2 represents the length of the sequence (like omega_2)
W2 = 50

# --- Step 1: Construct a family of functions F ---
# F will be a list of W2 functions. Each function is a list of W1 integers.
# We will construct F to be "increasing" and to diagonalize against
# potential bounding functions.

F = []

# Create the first function, f_0
f0 = [0] * W1
F.append(f0)

# Create the rest of the W2-1 functions by induction
for alpha in range(1, W2):
    f_prev = F[alpha-1]
    f_curr = [0] * W1

    # Make the new function strictly greater than the previous one everywhere
    for gamma in range(W1):
        f_curr[gamma] = f_prev[gamma] + 1

    # This is the diagonalization step. We make f_curr grow exceptionally
    # fast at a specific coordinate to "ruin" a potential bound.
    # A full proof would diagonalize against all possible bounds.
    # Here, we just add extra growth at a cycling coordinate.
    diag_coord = alpha % W1
    f_curr[diag_coord] += alpha # Add extra growth
    
    F.append(f_curr)

# --- Step 2: Define a large subset and a potential bounding function ---

# Let's take a "large" (uncountable) subset of the indices of F.
# Here, we'll take every other function starting from index 20.
X = range(20, W2, 2)

# Let's define an arbitrary function g that we'll test as a potential bound.
# This g grows linearly, which one might think could bound the family.
g = [5 * i for i in range(W1)]

print(f"Testing if the subset of functions F[beta] for beta in {list(X)} can be bounded by g.")
print(f"g = {g}\n")


# --- Step 3: Show that g fails to bound the subset X ---

bound_found = True
for beta in X:
    f_beta = F[beta]
    # Check if g bounds f_beta (i.e., g[gamma] > f_beta[gamma] for all gamma)
    is_bounded = True
    for gamma in range(W1):
        if f_beta[gamma] >= g[gamma]:
            # Found a coordinate where g fails to bound f_beta
            print(f"Found a counterexample with f_{beta}:")
            print(f"f_{beta} = {f_beta}")
            print(f"At coordinate gamma = {gamma}:")
            
            # Output the required "equation" format
            print(f"f_{beta}({gamma}) = {f_beta[gamma]}")
            print(f"g({gamma}) = {g[gamma]}")
            print(f"Result: f_{beta}({gamma}) >= g({gamma}), so g fails to bound f_{beta} and thus fails to bound the subset.")
            
            is_bounded = False
            bound_found = False
            break # Move to the next beta, or just stop
    if not is_bounded:
        break # We've already found a failure, so we can stop searching

if bound_found:
    print("Could not find a counterexample. The function g appears to bound the subset.")
