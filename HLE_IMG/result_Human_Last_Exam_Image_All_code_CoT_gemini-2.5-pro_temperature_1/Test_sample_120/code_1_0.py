# Data estimated from the plots at T = 375 K
# For N1 (non-methylated), the relaxation time of ring 1 (blue diamond)
tau_N1_ring1_at_375K = 8  # ns

# For M1 (methylated), the relaxation time of ring 1 (blue diamond)
tau_M1_ring1_at_375K = 15 # ns

print("Analysis of Relaxation Time at T = 375 K:")
print(f"Relaxation time for Ring 1 in N1 (non-methylated): {tau_N1_ring1_at_375K} ns")
print(f"Relaxation time for Ring 1 in M1 (methylated):    {tau_M1_ring1_at_375K} ns")

# Compare the values
if tau_M1_ring1_at_375K > tau_N1_ring1_at_375K:
    print("\nConclusion from data:")
    print("The relaxation time for the methylated ring is longer than for the non-methylated ring.")
    print("This means the addition of the methyl group increases the relaxation time, indicating slower rotational dynamics.")
else:
    print("\nConclusion from data:")
    print("The relaxation time for the methylated ring is not longer than for the non-methylated ring.")
