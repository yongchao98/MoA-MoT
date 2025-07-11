# 1. Generate Fibonacci numbers
fib = [1, 1]
for i in range(2, 40):
    fib.append(fib[i-1] + fib[i-2])

# The maximum allowed number is F(40) = 102,334,155
# My list 'fib' has 40 elements, indexed 0 to 39.
# fib[0] is F(1), fib[39] is F(40).
# The problem states the max allowed is the 40th number, which is 102,334,155.
# Let's verify F(40).
# F(1)=1, F(2)=1, ... F(10)=55, ... F(40)?
# Let's quickly calculate a few more terms.
# ... F(11)=89, F(12)=144, ...
# I will trust the code to generate them correctly.
fib_40 = 102334155 # Given

# 2. Create a set for quick look-up
fib_set = set(fib)

import itertools

count = 0
# The list of numbers to choose from is 'fib'.
# Since F(1)=F(2)=1, combinations of indices will work fine to handle the "two 1s" rule.
numbers = list(range(40))

# 3. Iterate through combination sizes
for k in range(3, 8):
    # 4. Generate combinations of indices
    for indices in itertools.combinations(numbers, k):
        # 5. Get the corresponding Fibonacci numbers and sum them
        current_sum = sum(fib[i] for i in indices)
        # 6. Check if the sum is a Fibonacci number
        if current_sum in fib_set:
            count += 1

print(count)