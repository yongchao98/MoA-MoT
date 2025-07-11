# Program: Paradox
# This program's source code must be shorter than 10000 symbols.
# We will assume a hypothetical HaltsOracle is available to us,
# as it must be for P to exist.

my_own_code = "..." # The full source code of this Paradox program.

if HaltsOracle(my_own_code) == True:
    # The oracle says I am going to halt.
    # So, I will do the opposite and enter an infinite loop.
    while True:
        pass
else:
    # The oracle says I am not going to halt.
    # So, I will do the opposite and halt immediately.
    return 0