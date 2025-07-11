import math

# --- Given Parameters ---
# For fast convolution methods like Overlap-Add/Save, the length of the
# shorter sequence (the filter, L) must be less than the DFT size (N).
# Therefore, the length-90 sequence is L and the length-1200 sequence is M.
L = 90
M = 1200
N = 128

print(f"Given parameters:")
print(f"Shorter sequence length (L): {L}")
print(f"Longer sequence length (M): {M}")
print(f"DFT/IDFT size (N): {N}\n")

# --- Overlap-Add Method ---
# The number of (DFT+IDFT) operations is determined by the number of blocks
# the longer sequence is broken into.

print("--- Overlap-Add Implementation ---")
# 1. Determine the block size for the long sequence (M).
# To avoid time-domain aliasing from the circular convolution, the block size
# (L_section) must satisfy: L_section + L - 1 <= N.
# We choose the maximum possible block size for efficiency.
# Formula: Block Size = N - L + 1
L_section_add = N - L + 1
print(f"1. Determine the block size for the long sequence.")
print(f"   Formula: Block Size = N - L + 1")
print(f"   Calculation: {N} - {L} + 1 = {L_section_add}")

# 2. Determine the number of blocks.
# The long sequence is segmented into non-overlapping blocks of size L_section_add.
# Formula: Num Blocks = ceil(M / Block Size)
num_blocks_add = math.ceil(M / L_section_add)
print(f"2. Determine the number of blocks.")
print(f"   Formula: Num Blocks = ceil(M / Block Size)")
print(f"   Calculation: ceil({M} / {L_section_add}) = {num_blocks_add}")
print(f"   Each block requires one DFT and one IDFT. So, {num_blocks_add} (DFT+IDFT) operations are needed.\n")


# --- Overlap-Save Method ---
print("--- Overlap-Save Implementation ---")
# 1. Determine the number of new samples processed per block (step size).
# Blocks are of size N and overlap by L-1 samples. The number of new samples
# in each block is N - (L-1).
# Formula: Step Size = N - L + 1
step_size_save = N - L + 1
print(f"1. Determine the number of new samples per block (step size).")
print(f"   Formula: Step Size = N - L + 1")
print(f"   Calculation: {N} - {L} + 1 = {step_size_save}")

# 2. Determine the number of blocks.
# The first block uses N samples. Each subsequent block advances by 'step_size_save'.
# We need to cover all M samples.
# Num Blocks = 1 (for first block) + ceil((M - N) / Step Size)
if M <= N:
    num_blocks_save = 1
else:
    remaining_samples = M - N
    num_subsequent_blocks = math.ceil(remaining_samples / step_size_save)
    num_blocks_save = 1 + num_subsequent_blocks
print(f"2. Determine the number of blocks.")
print(f"   The first block covers {N} samples. Remaining samples = {M} - {N} = {M - N}.")
print(f"   Each subsequent block covers {step_size_save} new samples.")
print(f"   Num subsequent blocks = ceil({M - N} / {step_size_save}) = {math.ceil((M - N) / step_size_save)}")
print(f"   Total Num Blocks = 1 + {math.ceil((M - N) / step_size_save)} = {num_blocks_save}")
print(f"   Each block requires one DFT and one IDFT. So, {num_blocks_save} (DFT+IDFT) operations are needed.\n")

# --- Final Answer ---
print("--- Summary ---")
print(f"Overlap-Add requires: {num_blocks_add} (DFT+IDFT) operations.")
print(f"Overlap-Save requires: {num_blocks_save} (DFT+IDFT) operations.")