# --- Parameters ---
# Approximate length of a TCR beta chain mRNA, excluding the poly(A) tail.
TCR_TRANSCRIPT_LENGTH = 1300  # in bases

# Approximate start position of the CDR3 region, measured from the 5' end.
CDR3_START_FROM_5_PRIME = 350 # in bases

# The length of the sequencing read used to sequence the transcript body.
READ_2_LENGTH = 225 # in bases

# Approximate start of the TCR Constant (C) region, measured from the 3' end.
# This is our new target for the specific capture oligo.
C_REGION_PRIMING_SITE_FROM_3_PRIME = 400 # in bases

# --- Calculations ---

# 1. Calculate the position of the CDR3 region from the 3' end, which is where sequencing starts.
cdr3_start_from_3_prime = TCR_TRANSCRIPT_LENGTH - CDR3_START_FROM_5_PRIME

# --- Scenario 1: Original method with Poly(dT) priming ---
print("--- Scenario 1: Original Poly(dT) Priming ---")
# The read starts at the 3' end (position 0) and extends for READ_2_LENGTH.
poly_dt_read_coverage_end = READ_2_LENGTH
print(f"The read sequences from base 0 to base {poly_dt_read_coverage_end} from the 3' end.")
print(f"The CDR3 region begins at base {cdr3_start_from_3_prime} from the 3' end.")

# Check if the CDR3 is covered by the read.
if poly_dt_read_coverage_end >= cdr3_start_from_3_prime:
    print("Result: Success! The CDR3 region is sequenced.\n")
else:
    print(f"Result: Failure. The read ends at base {poly_dt_read_coverage_end}, but the CDR3 region starts at base {cdr3_start_from_3_prime}. The read is too short by {cdr3_start_from_3_prime - poly_dt_read_coverage_end} bases.\n")


# --- Scenario 2: Proposed method with TCR-specific priming (Option A) ---
print("--- Scenario 2: Proposed TCR-Specific Priming in C-Region ---")
# The read now starts at the new priming site in the constant region.
specific_read_coverage_start = C_REGION_PRIMING_SITE_FROM_3_PRIME
specific_read_coverage_end = C_REGION_PRIMING_SITE_FROM_3_PRIME + READ_2_LENGTH
print(f"The new priming site is at base {specific_read_coverage_start} from the 3' end.")
print(f"The read now sequences the region from base {specific_read_coverage_start} to base {specific_read_coverage_end} from the 3' end.")
print(f"The CDR3 region begins at base {cdr3_start_from_3_prime} from the 3' end.")

# Check if the CDR3 is covered by the new read position.
if specific_read_coverage_end >= cdr3_start_from_3_prime and specific_read_coverage_start <= cdr3_start_from_3_prime:
    print("Result: Success! The CDR3 region is now within the sequenced portion of the transcript.\n")
else:
    print("Result: Failure. The CDR3 region is not covered. (Check priming site selection).\n")
