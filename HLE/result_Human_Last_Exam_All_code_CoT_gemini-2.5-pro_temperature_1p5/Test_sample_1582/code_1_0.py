import random

def simulate_ssrw_excursions(num_excursions):
    """
    Simulates a Simple Symmetric Random Walk on Z and records
    the duration of excursions from state 0.
    
    Args:
        num_excursions (int): The number of excursions to simulate.
    """
    
    print("This simulation demonstrates a non-positive recurrent chain satisfying the problem's conditions.")
    print("Chain: Simple Symmetric Random Walk on Z. f(x) = |x|. A = {0}.")
    print("The condition is E[f(X_1)|X_0=x] - f(x) >= 0 for x != 0.")
    print("Let's check for x = 5:")
    x = 5
    f_x_plus_1 = abs(x + 1)
    f_x_minus_1 = abs(x - 1)
    f_x = abs(x)
    drift = 0.5 * f_x_plus_1 + 0.5 * f_x_minus_1 - f_x
    print(f"0.5 * f({x+1}) + 0.5 * f({x-1}) - f({x}) = 0.5 * {f_x_plus_1} + 0.5 * {f_x_minus_1} - {f_x} = {drift:.1f}")
    print("-" * 30)

    current_pos = 0
    excursions_completed = 0
    excursion_lengths = []
    
    print(f"Simulating {num_excursions} excursions from 0...")
    
    while excursions_completed < num_excursions:
        # Start an excursion by moving from 0
        if random.random() < 0.5:
            current_pos = 1
        else:
            current_pos = -1
        
        current_excursion_length = 1
        
        # Continue the walk until it returns to 0
        while current_pos != 0:
            if random.random() < 0.5:
                current_pos += 1
            else:
                current_pos -= 1
            current_excursion_length += 1
            
        # Excursion has ended
        excursions_completed += 1
        excursion_lengths.append(current_excursion_length)
        
    print(f"Finished simulation.")
    print(f"Durations of the first {min(20, num_excursions)} excursions: {excursion_lengths[:20]}")
    print(f"Average excursion length: {sum(excursion_lengths) / len(excursion_lengths):.2f}")
    print(f"Maximum excursion length observed: {max(excursion_lengths)}")
    print("\nNote: For a null recurrent chain, there is no finite average excursion time.")
    print("If you run this simulation multiple times or for more excursions, you'll see very large maximums appear,")
    print("and the average will not converge to a stable value, indicating non-positive recurrence.")

if __name__ == '__main__':
    # Running with a decent number of excursions to see the effect
    simulate_ssrw_excursions(num_excursions=1000)