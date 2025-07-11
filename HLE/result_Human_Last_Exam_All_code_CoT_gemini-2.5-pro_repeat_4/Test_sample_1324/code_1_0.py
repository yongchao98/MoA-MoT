import math

def simulate_chocolates(initial_chocolates, steps):
    """
    Simulates the chocolate passing game and prints the evolution of key metrics.

    Args:
        initial_chocolates (list): The initial number of chocolates for each person.
        steps (int): The maximum number of minutes to simulate.
    """
    chocolates = list(initial_chocolates)
    n = len(chocolates)

    print(f"Let's simulate with n = {n} people.")
    print(f"The initial distribution of chocolates is c^0 = {chocolates}.")
    print("-" * 30)

    for i in range(steps + 1):
        h = max(chocolates)
        l = min(chocolates)
        d = h - l

        print(f"Minute i={i}:")
        print(f"  Chocolates c^{i}: {chocolates}")
        print(f"  Max h^{i} = {h}, Min l^{i} = {l}, Difference d^{i} = {d}")
        
        if d == 0:
            print("\nThe system has reached a stable state.")
            print("Let's check Statements 1 and 3 for this stable step i={}:".format(i))
            print("Statement 1 requires d^{}>d^{}, which is 0 < 0. This is false.".format(i+1, i))
            print("Statement 3 requires l^{}>l^{}, which is {} > {}. This is false.".format(i+1, i, l, l))
            break
        
        print("-" * 30)

        new_chocolates = [0] * n
        for k in range(n):
            # Person k receives from person k-1 (circularly)
            prev_person_idx = (k - 1 + n) % n
            
            c_sum = chocolates[k] + chocolates[prev_person_idx]
            half_sum = c_sum / 2
            
            # The new amount is the smallest even integer >= half_sum.
            # This is equivalent to 2 * ceil(sum / 4)
            new_amount = math.ceil(half_sum)
            if new_amount % 2 != 0:
                new_amount += 1
            new_chocolates[k] = new_amount

        chocolates = new_chocolates

# Example run
n_people = 4
initial_dist = [2, 4, 6, 12] 
simulation_steps = 10

simulate_chocolates(initial_dist, simulation_steps)