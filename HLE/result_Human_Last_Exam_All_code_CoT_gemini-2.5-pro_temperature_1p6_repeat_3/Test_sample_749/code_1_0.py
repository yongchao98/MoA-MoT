import math

# Step 1: Define jump probabilities in the limit h->0.
# As h->0, all sites become blue with probability 1.
# For a blue site, P(left) = 1/5, P(right) = 4/5.
p_L = 1/5
p_R = 4/5

print(f"In the h->0 limit, the walk is biased with:")
print(f"Probability of jumping left, p_L = {p_L}")
print(f"Probability of jumping right, p_R = {p_R}")
print("-" * 30)

# Step 2: Calculate hitting probabilities for this biased random walk.
# The walk has a drift to the right since p_R > p_L.
# The probability of reaching 0 from any site k < 0 is 1.
pi_minus_1 = 1.0

# The probability of reaching 0 from site k=1 is p_L / p_R.
pi_1 = p_L / p_R

print(f"The walk is transient to the right.")
print(f"The probability of reaching site 0 from site -1 is pi_(-1) = {pi_minus_1}")
print(f"The probability of reaching site 0 from site 1 is pi_1 = p_L / p_R = {p_L/p_R:.2f}")
print("-" * 30)


# Step 3: Calculate m_0, the expected number of returns to 0 in the h->0 limit.
# A particle at 0 jumps to -1 (prob p_L) or +1 (prob p_R).
# m_0 is the probability of returning from that new position.
m_0 = p_L * pi_minus_1 + p_R * pi_1

print(f"The mean number of returning descendants in the h->0 limit is m_0.")
print(f"This is the probability for a particle starting at 0 to ever return.")
print(f"m_0 = p_L * pi_(-1) + p_R * pi_1")
print(f"m_0 = {p_L} * {pi_minus_1} + {p_R} * {pi_1:.2f}")
print(f"m_0 = {p_L * pi_minus_1} + {p_R * pi_1}")
print(f"m_0 = {m_0}")
print("-" * 30)

# Step 4: Conclude based on the value of m_0.
if m_0 < 1:
    print(f"Since m_0 = {m_0} is less than 1, the branching process of returns to site 0 is subcritical in the h->0 limit.")
    print("This implies that for any sufficiently small h > 0, the process is subcritical.")
    print("A subcritical branching process dies out with probability 1.")
    print("Therefore, the number of particles visiting site 0 is finite with probability 1.")
    final_limit = 0
else:
    # This case won't be reached with these numbers.
    final_limit = "Calculation indicates a non-zero limit."

print(f"\nThe final result is:")
print(f"lim_{{h->0}} P[site 0 is visited by infinitely many different particles] = {final_limit}")