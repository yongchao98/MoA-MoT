import math

def get_divisors(n):
    """Returns a list of all divisors of a number n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return list(divs)

def solve_subset_sum(numbers, target_sum):
    """
    Solves the subset sum problem using dynamic programming.
    Returns True if a subset of 'numbers' sums to 'target_sum', False otherwise.
    """
    dp = {0}
    for num in numbers:
        dp.update({x + num for x in dp})
        if target_sum in dp:
            return True
    return target_sum in dp

def is_zumkeller(n):
    """
    Checks if a number n is a Zumkeller number.
    A number n is a Zumkeller number if its divisors can be partitioned
    into two sets with equal sums.
    """
    if n == 0:
        return False
        
    divisors = get_divisors(n)
    sigma_n = sum(divisors)
    
    # The sum of divisors must be even to be split into two equal integer sums.
    if sigma_n % 2 != 0:
        return False
        
    target = sigma_n // 2
    
    return solve_subset_sum(divisors, target)

def find_next_number_in_sequence():
    """
    Finds and prints the next number in the sequence.
    """
    last_number = 2352
    
    print("The pattern in the sequence is that each number is a 'Zumkeller number'.")
    print("A Zumkeller number is a number whose divisors can be partitioned into two sets with equal sums.")
    print(f"The last number given is {last_number}. We will search for the next Zumkeller number starting from {last_number + 1}.")
    
    current_number = last_number + 1
    while True:
        if is_zumkeller(current_number):
            next_number = current_number
            break
        current_number += 1
        
    difference = next_number - last_number
    
    print(f"\nThe next Zumkeller number after {last_number} is {next_number}.")
    print("The final equation is:")
    print(f"{last_number} + {difference} = {next_number}")

if __name__ == "__main__":
    find_next_number_in_sequence()
    # The calculated next number is 2360.
    # To follow the output format, we will output this number within <<<>>>.
    # print("<<<2360>>>") is commented out to avoid printing it twice in execution.