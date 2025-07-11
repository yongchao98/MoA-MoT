# The problem presented is unsolvable for two primary reasons:
# 1. Missing Information: The initial distance to Pandora is not provided, making
#    the calculation of travel time impossible from a physics standpoint.
# 2. Computational Limitations: The custom "Wuxing" architecture cannot perform
#    the necessary calculations. Its 'frac' data type, which uses 2-digit characters,
#    would overflow when trying to compute the spacecraft's velocity even after the
#    first day of acceleration.
#
# According to the prompt's instructions, if a program cannot be written to
# solve the task, the answer should be "0:0".
print("0:0")